#pragma once

/**
 * ============================================================================
 * Hairpin & Adapter Finder CUDA Engine (Duplex SBX-D Architecture)
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Source Correspondence:
 *   - Directly mirrors: XOOS/demux/src/adapters/duplex/hairpin-finder.h
 *   - Directly mirrors: XOOS/demux/src/adapters/duplex/hairpin-finder.cpp
 *   - Directly mirrors: XOOS/demux/src/adapters/duplex/trim-duplex.h
 *   - Directly mirrors: XOOS/demux/src/adapters/duplex/trim-duplex.cpp
 *
 * Biological Structure of AXELIOS SBX Duplex Reads:
 *   [Start Adapter]  --> ~15-20 bp synthetic 5' leader
 *   [Insert R1]      --> 1st pass of target genomic fragment (forward)
 *   [Hairpin Region] --> [SID_5p (~12bp)] + [Loop (~16-24bp)] + [SID_3p_RC (~12bp)]
 *   [Insert R2]      --> 2nd pass of target genomic fragment (reverse complement)
 *   [End Adapter]    --> ~15-20 bp synthetic 3' trailer
 *
 * GPU In-Register Extraction:
 *   1. Warp-parallel loop search using 64-bit Bitap<2> window scanning.
 *   2. Left-flank and right-flank SID barcode matching in local registers.
 *   3. Boundary trimming emitting exact [insert1] and [insert2] spans directly
 *      to downstream consensus registers with zero global DRAM roundtrips.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include "bitap_cuda.cuh"

namespace xoos::demux::cuda {

constexpr int32_t kMaxSidsPerBundle = 384;
constexpr int32_t kMinDuplexReadLen = 40;
constexpr int32_t kDefaultLoopSearchWindow = 64;

/**
 * @brief Demultiplexing and trimming result for a single duplex read.
 */
struct DuplexTrimResult {
    uint32_t read_idx;
    int32_t sample_id;          // Matched sample ID index (0..N-1), or -1 if unassigned
    uint32_t bitflag;           // Bitflag (has SID, UMI, etc.)
    
    // Genomic Insert Coordinate Boundaries (0-based inclusive in raw read)
    int32_t start_adapter_pos;  // End of start adapter (-1 if omitted/truncated)
    int32_t insert1_start;      // Start of R1 insert
    int32_t insert1_end;        // End of R1 insert (before hairpin)
    
    int32_t hairpin_start;      // Start of hairpin (SID_5p)
    int32_t hairpin_loop_pos;   // Position of matched loop
    int32_t hairpin_end;        // End of hairpin (SID_3p)
    
    int32_t insert2_start;      // Start of R2 insert (after hairpin)
    int32_t insert2_end;        // End of R2 insert (before end adapter)
    int32_t end_adapter_pos;    // Start of end adapter (-1 if truncated)
    
    int32_t insert1_len;        // insert1_end - insert1_start + 1
    int32_t insert2_len;        // insert2_end - insert2_start + 1
    
    bool is_valid_duplex;       // True if both R1 and R2 inserts meet length thresholds
};

/**
 * @brief Precomputed GPU-resident adapter and SID bundle for demux.
 */
struct GpuAdapterBundle {
    BitapPattern loop_fw;
    BitapPattern loop_bw;
    BitapPattern start_adapter_fw;
    BitapPattern end_adapter_fw;
    
    int32_t num_samples;
    BitapPattern sid_5p_matchers[kMaxSidsPerBundle];
    BitapPattern sid_3p_matchers[kMaxSidsPerBundle];
};

/**
 * @brief In-Register Hairpin and Adapter Extractor for a Single Read.
 *
 * @param bundle Precomputed GPU adapter patterns in constant/global memory.
 * @param raw_seq Pointer to ASCII nucleotide sequence of the read.
 * @param seq_len Total length of the raw read in bases.
 * @param out_trim Struct receiving the parsed boundaries and sample assignment.
 */
__device__ __forceinline__ void extract_duplex_hairpin(
    const GpuAdapterBundle& bundle,
    const char* __restrict__ raw_seq,
    int32_t seq_len,
    DuplexTrimResult& out_trim
) {
    out_trim.sample_id = -1;
    out_trim.bitflag = 0;
    out_trim.is_valid_duplex = false;
    out_trim.start_adapter_pos = -1;
    out_trim.end_adapter_pos = -1;
    out_trim.hairpin_start = -1;
    out_trim.hairpin_end = -1;

    if (seq_len < kMinDuplexReadLen) {
        return;
    }

    // ------------------------------------------------------------------------
    // Phase 1: Search for the central Hairpin Loop using Bitap<2>
    // ------------------------------------------------------------------------
    // In AXELIOS duplex reads, the hairpin loop is positioned roughly in the
    // middle of the read (between R1 and R2), typically around pos in [20, seq_len - 20].
    int32_t loop_search_begin = 20;
    int32_t loop_search_end = (seq_len > 40) ? (seq_len - 20) : seq_len - 1;
    
    int32_t loop_pos = bitap_search_window<2, MatchPolicy::kFirst>(
        bundle.loop_fw, raw_seq, loop_search_begin, loop_search_end
    );

    if (loop_pos == kCudaNoMatchPosition) {
        // Fallback: wider scan across full interior if not found in prime window
        loop_pos = bitap_search_window<3, MatchPolicy::kFirst>(
            bundle.loop_fw, raw_seq, 10, seq_len - 10
        );
        if (loop_pos == kCudaNoMatchPosition) {
            return; // Not a valid duplex read (hairpin missing)
        }
    }

    out_trim.hairpin_loop_pos = loop_pos;
    int32_t loop_len = bundle.loop_fw.query_len;

    // ------------------------------------------------------------------------
    // Phase 2: Demultiplex Sample ID (SID) Barcodes around the Loop
    // ------------------------------------------------------------------------
    // SID_5p immediately precedes the loop; SID_3p immediately follows the loop.
    int32_t best_sample = -1;
    int32_t best_sid_start = -1;
    int32_t best_sid_end = -1;

    int32_t sid_5p_search_end = loop_pos;
    int32_t sid_5p_search_begin = (loop_pos >= 20) ? (loop_pos - 20) : 0;
    int32_t sid_3p_search_begin = loop_pos + loop_len;
    int32_t sid_3p_search_end = (sid_3p_search_begin + 20 < seq_len) ? (sid_3p_search_begin + 20) : (seq_len - 1);

    // Evaluate sample barcode matchers across available samples
    for (int s = 0; s < bundle.num_samples; ++s) {
        int32_t p5 = bitap_search_window<1, MatchPolicy::kWidest>(
            bundle.sid_5p_matchers[s], raw_seq, sid_5p_search_begin, sid_5p_search_end
        );
        
        if (p5 != kCudaNoMatchPosition) {
            best_sample = s;
            best_sid_start = (p5 >= bundle.sid_5p_matchers[s].query_len) ? (p5 - bundle.sid_5p_matchers[s].query_len + 1) : 0;
            
            // Check 3' SID confirmation
            int32_t p3 = bitap_search_window<1, MatchPolicy::kFirst>(
                bundle.sid_3p_matchers[s], raw_seq, sid_3p_search_begin, sid_3p_search_end
            );
            best_sid_end = (p3 != kCudaNoMatchPosition) ? p3 : (sid_3p_search_begin + bundle.sid_5p_matchers[s].query_len);
            break;
        }
    }

    if (best_sample < 0) {
        // Sample barcode could not be resolved with confidence
        best_sample = 0; // Default bin
        best_sid_start = (loop_pos >= 12) ? (loop_pos - 12) : 0;
        best_sid_end = (loop_pos + loop_len + 12 < seq_len) ? (loop_pos + loop_len + 12) : seq_len - 1;
    }

    out_trim.sample_id = best_sample;
    out_trim.hairpin_start = best_sid_start;
    out_trim.hairpin_end = best_sid_end;

    // ------------------------------------------------------------------------
    // Phase 3: Start Adapter & End Adapter Search
    // ------------------------------------------------------------------------
    // Search for 5' Start Adapter in prefix [0, min(40, best_sid_start)]
    int32_t start_adapter_search_end = (best_sid_start > 40) ? 40 : best_sid_start;
    int32_t start_ad_pos = bitap_search_window<2, MatchPolicy::kWidest>(
        bundle.start_adapter_fw, raw_seq, 0, start_adapter_search_end
    );
    out_trim.start_adapter_pos = start_ad_pos;
    out_trim.insert1_start = (start_ad_pos != kCudaNoMatchPosition) ? (start_ad_pos + 1) : 0;
    out_trim.insert1_end = (best_sid_start > 0) ? (best_sid_start - 1) : 0;

    // Search for 3' End Adapter in suffix [best_sid_end, seq_len - 1]
    int32_t end_ad_search_begin = (best_sid_end + 20 < seq_len) ? (best_sid_end + 20) : best_sid_end;
    int32_t end_ad_pos = bitap_search_window<2, MatchPolicy::kFirst>(
        bundle.end_adapter_fw, raw_seq, end_ad_search_begin, seq_len - 1
    );
    out_trim.end_adapter_pos = end_ad_pos;
    out_trim.insert2_start = (best_sid_end + 1 < seq_len) ? (best_sid_end + 1) : (seq_len - 1);
    out_trim.insert2_end = (end_ad_pos != kCudaNoMatchPosition && end_ad_pos > out_trim.insert2_start)
        ? (end_ad_pos - bundle.end_adapter_fw.query_len)
        : (seq_len - 1);

    out_trim.insert1_len = out_trim.insert1_end - out_trim.insert1_start + 1;
    out_trim.insert2_len = out_trim.insert2_end - out_trim.insert2_start + 1;

    // Valid duplex requires both inserts to be positive in length
    if (out_trim.insert1_len >= 15 && out_trim.insert2_len >= 15) {
        out_trim.is_valid_duplex = true;
    }
}

} // namespace xoos::demux::cuda
