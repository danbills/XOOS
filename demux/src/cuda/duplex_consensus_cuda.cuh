#pragma once

/**
 * ============================================================================
 * Duplex Consensus Collapser & Quality Solver CUDA Engine
 * ============================================================================
 *
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Source Correspondence:
 *   - Directly mirrors: XOOS/demux/src/utility/alignment-util.h
 *   - Directly mirrors: XOOS/demux/src/utility/alignment-util.cpp
 *   - Directly mirrors: XOOS/demux/src/task/alignment.h
 *   - Directly mirrors: XOOS/demux/src/task/alignment.cpp
 *   - Directly mirrors: XOOS/demux/src/task/formatter.cpp
 *
 * Consensus Alignment & Quality Recalibration Mechanics:
 *   In AXELIOS SBX Duplex reads, R1 and R2 represent the two complementary
 *   strands of the same double-stranded DNA molecule:
 *     R1: forward strand sequence [insert1_start .. insert1_end]
 *     R2: reverse strand sequence [insert2_start .. insert2_end]
 *
 *   Pairwise comparison aligns R1[i] with the reverse-complement of R2:
 *     R2_RC[i] = complement(R2[insert2_end - i])
 *
 *   1. Concordant Base (R1[i] == R2_RC[i]):
 *      Both strands agree. Base error probability drops multiplicatively:
 *        Q_consensus = min(60, Q1 + Q2)  --> Yields Q40 to Q60 ultra-high confidence.
 *   2. Discordant Base (R1[i] != R2_RC[i]):
 *      Strands disagree due to single-strand PCR or sequencing error:
 *        Q_consensus = max(0, |Q1 - Q2|)
 *      Emits IUPAC ambiguity code for the YC consensus tag (e.g. A+G -> R, C+T -> Y).
 *   3. Epigenetic Methylation Conversion (CpG Dinucleotides):
 *      Scans adjacent CG sites to classify 5mC methylation states:
 *        'Z' = Fully methylated (both strands)
 *        'U' = Hemi-methylated on R1
 *        'u' = Hemi-methylated on R2
 *        'z' = Unmethylated CpG
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <cuda_runtime.h>
#include <cstdint>
#include "hairpin_cuda.cuh"

namespace xoos::demux::cuda {

constexpr int32_t kMaxConsensusLen = 128;

/**
 * @brief Complement of DNA character (A<->T, C<->G).
 */
__device__ __forceinline__ char complement_dna(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 'T';
        case 'c': return 'G';
        case 'g': return 'C';
        case 't': return 'A';
        default:  return 'N';
    }
}

/**
 * @brief IUPAC Ambiguity Code calculation for discordant base pairs.
 * Mirrors IUPAC encoding in XOOS alignment-util.
 */
__device__ __forceinline__ char get_iupac_code(char b1, char b2) {
    char u1 = (b1 >= 'a' && b1 <= 'z') ? (b1 - 32) : b1;
    char u2 = (b2 >= 'a' && b2 <= 'z') ? (b2 - 32) : b2;

    if (u1 == u2) return u1;
    if ((u1 == 'A' && u2 == 'G') || (u1 == 'G' && u2 == 'A')) return 'R'; // PuRine
    if ((u1 == 'C' && u2 == 'T') || (u1 == 'T' && u2 == 'C')) return 'Y'; // PYrimidine
    if ((u1 == 'A' && u2 == 'C') || (u1 == 'C' && u2 == 'A')) return 'M'; // aMino
    if ((u1 == 'G' && u2 == 'T') || (u1 == 'T' && u2 == 'G')) return 'K'; // Keto
    if ((u1 == 'C' && u2 == 'G') || (u1 == 'G' && u2 == 'C')) return 'S'; // Strong
    if ((u1 == 'A' && u2 == 'T') || (u1 == 'T' && u2 == 'A')) return 'W'; // Weak
    return 'N';
}

/**
 * @brief Formatted consensus read record ready for FASTQ emission.
 */
struct ConsensusReadResult {
    int32_t sample_id;
    int32_t consensus_len;
    char consensus_seq[kMaxConsensusLen];
    char consensus_qual[kMaxConsensusLen];
    char yc_tag[kMaxConsensusLen];
    char xm_tag[kMaxConsensusLen];
    
    int32_t concordant_bases;
    int32_t discordant_bases;
    float concordance_rate;
};

/**
 * @brief Collapse Duplex R1 and R2 into Intramolecular Consensus in Registers.
 *
 * @param raw_seq Pointer to raw nucleotide sequence of the read.
 * @param raw_qual Pointer to raw Phred quality string of the read.
 * @param trim Trimming boundaries computed by extract_duplex_hairpin.
 * @param out_res Struct receiving the consensus sequence, quality, and tags.
 */
__device__ __forceinline__ void collapse_duplex_consensus(
    const char* __restrict__ raw_seq,
    const char* __restrict__ raw_qual,
    const DuplexTrimResult& trim,
    ConsensusReadResult& out_res
) {
    out_res.sample_id = trim.sample_id;
    out_res.consensus_len = 0;
    out_res.concordant_bases = 0;
    out_res.discordant_bases = 0;
    out_res.concordance_rate = 0.0f;

    if (!trim.is_valid_duplex) {
        return;
    }

    const int32_t r1_start = trim.insert1_start;
    const int32_t r1_len = trim.insert1_len;
    const int32_t r2_end = trim.insert2_end;
    const int32_t r2_len = trim.insert2_len;

    // Consensus length is the minimum overlapping insert length
    const int32_t cons_len = (r1_len < r2_len) ? r1_len : r2_len;
    if (cons_len <= 0 || cons_len >= kMaxConsensusLen) {
        return;
    }

    out_res.consensus_len = cons_len;

    // Phase 1: Pairwise base comparison and quality score recalibration
    for (int32_t i = 0; i < cons_len; ++i) {
        char b1 = raw_seq[r1_start + i];
        int8_t q1 = static_cast<int8_t>(raw_qual[r1_start + i] - 33);
        if (q1 < 0) q1 = 0;

        // R2 base taken in reverse complement order
        char b2_raw = raw_seq[r2_end - i];
        char b2_rc = complement_dna(b2_raw);
        int8_t q2 = static_cast<int8_t>(raw_qual[r2_end - i] - 33);
        if (q2 < 0) q2 = 0;

        char final_base;
        int32_t final_q;
        char iupac = get_iupac_code(b1, b2_rc);

        if (b1 == b2_rc) {
            // Concordant base call
            final_base = b1;
            int32_t sum_q = q1 + q2;
            final_q = (sum_q > 60) ? 60 : sum_q;
            out_res.concordant_bases++;
        } else {
            // Discordant base call
            final_base = (q1 >= q2) ? b1 : b2_rc;
            int32_t diff_q = (q1 > q2) ? (q1 - q2) : (q2 - q1);
            final_q = (diff_q < 2) ? 2 : diff_q;
            out_res.discordant_bases++;
        }

        out_res.consensus_seq[i] = final_base;
        out_res.consensus_qual[i] = static_cast<char>(final_q + 33);
        out_res.yc_tag[i] = iupac;
        out_res.xm_tag[i] = '.';
    }

    out_res.consensus_seq[cons_len] = '\0';
    out_res.consensus_qual[cons_len] = '\0';
    out_res.yc_tag[cons_len] = '\0';
    out_res.xm_tag[cons_len] = '\0';

    if (cons_len > 0) {
        out_res.concordance_rate = static_cast<float>(out_res.concordant_bases) / static_cast<float>(cons_len);
    }

    // Phase 2: CpG Methylation Status Detection
    // Scans 2 bases at a time for CpG patterns:
    //   'Z' = Fully methylated
    //   'U' = Hemi-methylated (R1)
    //   'u' = Hemi-methylated (R2)
    //   'z' = Unmethylated
    for (int32_t i = 0; i < cons_len - 1; ++i) {
        char c1 = out_res.consensus_seq[i];
        char c2 = out_res.consensus_seq[i + 1];
        char r2_1 = complement_dna(raw_seq[r2_end - i]);
        char r2_2 = complement_dna(raw_seq[r2_end - (i + 1)]);

        if (c1 == 'C' && c2 == 'G' && r2_1 == 'C' && r2_2 == 'G') {
            out_res.xm_tag[i] = 'z'; // Unmethylated CpG
        } else if (c1 == 'T' && c2 == 'G' && r2_1 == 'C' && r2_2 == 'A') {
            out_res.consensus_seq[i] = 'C';
            out_res.xm_tag[i] = 'Z'; // Fully methylated CpG
        } else if (c1 == 'T' && c2 == 'G' && r2_1 == 'C' && r2_2 == 'G') {
            out_res.consensus_seq[i] = 'C';
            out_res.xm_tag[i] = 'U'; // Hemi-methylated R1
        } else if (c1 == 'C' && c2 == 'G' && r2_1 == 'C' && r2_2 == 'A') {
            out_res.xm_tag[i] = 'u'; // Hemi-methylated R2
        }
    }
}

} // namespace xoos::demux::cuda
