#pragma once

/**
 * ============================================================================
 * Bitap CUDA Engine: 64-bit Register-Resident Approximate String Matcher
 * ============================================================================
 * 
 * Target Architecture: NVIDIA Blackwell (SM 120 / Compute Capability 12.0)
 *                      and CUDA 12.0+ / 13.x
 *
 * Source Correspondence:
 *   - Directly mirrors: XOOS/demux/src/sequence/matcher/bitap.h
 *   - Directly mirrors: XOOS/demux/src/sequence/matcher/bitap.cpp
 *   - Algorithm: Baeza-Yates–Gonnet (Shift-And / Bitap) with Myers-style
 *     dynamic programming bit-parallel level state transitions.
 *
 * GPU Architecture Optimizations:
 *   1. Zero DRAM Latency: State machine runs entirely in 64-bit ALU registers
 *      (uint64_t `current_state`, `previous_state`, `best_pos`).
 *   2. Constant/Register Alphabet: DNA alphabet masks (A, C, G, T) fit in
 *      4 x uint64_t registers per query (32 bytes), avoiding L1/L2 thrashing.
 *   3. Warp-Uniform Execution: Inner bitwise transitions (SHL, AND, OR)
 *      are branchless single-cycle instructions with zero thread divergence.
 *
 * Copyright (c) 2026 Roche Diagnostics / AXELIOS Open Source
 * ============================================================================
 */

#include <cuda_runtime.h>
#include <cstdint>
#include <cstddef>

namespace xoos::demux::cuda {

// Constant representing "not found" position (mirrors kNotFound in bitap.cpp)
constexpr int32_t kCudaNotFound = 0x7FFFFFFF;
constexpr int32_t kCudaNoMatchPosition = -1;
constexpr int32_t kCudaQueryWindowSize = 64;

enum class SearchDirection : uint8_t {
    kForward = 0,   // Returns right edge (e.g. 5' trim point after adapter)
    kBackward = 1   // Returns left edge  (e.g. 3' trim point before adapter)
};

enum class MatchPolicy : uint8_t {
    kFirst = 0,     // Early exit on first match at lowest edit distance
    kWidest = 1     // Scans full window, returning latest (widest) match
};

/**
 * @brief Base-to-index mapping for DNA characters (A=0, C=1, G=2, T=3, N=0)
 */
__host__ __device__ __forceinline__ uint8_t base_to_index(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 0; // Fallback / N
    }
}

/**
 * @brief 64-bit Bitap pattern mask descriptor.
 * Can reside in __constant__ memory or thread local registers.
 */
struct BitapPattern {
    uint64_t alphabet[4]; // [0]=A, [1]=C, [2]=G, [3]=T (inverted masks)
    uint64_t msb;         // (1ULL << (query_len - 1))
    int32_t query_len;
    SearchDirection direction;

    /**
     * @brief Host-side initialization of pattern masks.
     */
    __host__ __device__ void init(const char* query, int len, SearchDirection dir = SearchDirection::kForward) {
        query_len = len;
        direction = dir;
        msb = (len > 0 && len <= 64) ? (1ULL << (len - 1)) : 0ULL;

        for (int b = 0; b < 4; ++b) {
            alphabet[b] = 0ULL;
        }

        if (dir == SearchDirection::kBackward) {
            for (int i = 0; i < len; ++i) {
                char c = query[len - 1 - i];
                uint8_t idx = (c == 'A' || c == 'a') ? 0 :
                              (c == 'C' || c == 'c') ? 1 :
                              (c == 'G' || c == 'g') ? 2 :
                              (c == 'T' || c == 't') ? 3 : 0;
                alphabet[idx] |= (1ULL << i);
            }
        } else {
            for (int i = 0; i < len; ++i) {
                char c = query[i];
                uint8_t idx = (c == 'A' || c == 'a') ? 0 :
                              (c == 'C' || c == 'c') ? 1 :
                              (c == 'G' || c == 'g') ? 2 :
                              (c == 'T' || c == 't') ? 3 : 0;
                alphabet[idx] |= (1ULL << i);
            }
        }

        // Invert bit masks for bitap algorithm efficiency (~mask)
        for (int b = 0; b < 4; ++b) {
            alphabet[b] = ~alphabet[b];
        }
    }
};

/**
 * @brief Register-resident Bitap Approximate String Search on GPU.
 *
 * @tparam MaxDist Maximum edit distance allowed (0, 1, 2, 3, or 4).
 * @tparam Policy kFirst (earliest match) or kWidest (widest match).
 * @param pattern Precomputed query BitapPattern.
 * @param ref Pointer to reference DNA character buffer.
 * @param begin Starting offset in ref (inclusive).
 * @param end Ending offset in ref (inclusive).
 * @return 0-based match coordinate in ref, or kCudaNoMatchPosition if not found.
 */
template <uint32_t MaxDist, MatchPolicy Policy = MatchPolicy::kFirst>
__host__ __device__ __forceinline__ int32_t bitap_search_window(
    const BitapPattern& pattern,
    const char* __restrict__ ref,
    int32_t begin,
    int32_t end
) {
    const int32_t ref_length = 1 + end - begin;
    if (ref_length <= 0 || pattern.query_len <= 0) {
        return kCudaNoMatchPosition;
    }

    constexpr int32_t kNumLevels = MaxDist + 1;
    uint64_t current_state[kNumLevels];
    uint64_t previous_state[kNumLevels];
    int32_t best_position[kNumLevels];

    #pragma unroll
    for (int d = 0; d < kNumLevels; ++d) {
        current_state[d] = ~0ULL;
        previous_state[d] = ~0ULL;
        best_position[d] = kCudaNotFound;
    }

    // Process each character in the reference sequence
    for (int32_t i = 0; i < ref_length; ++i) {
        char base = ref[begin + i];
        uint8_t base_idx = base_to_index(base);
        uint64_t current_pattern_mask = pattern.alphabet[base_idx];

        // Process exact match (distance 0)
        previous_state[0] = current_state[0];
        current_state[0] = (current_state[0] << 1) | current_pattern_mask;

        if ((current_state[0] & pattern.msb) == 0) {
            if constexpr (Policy == MatchPolicy::kFirst) {
                // Early return on first exact match (zero edit distance)
                return (pattern.direction == SearchDirection::kBackward)
                    ? (end - i)
                    : (begin + i);
            } else {
                best_position[0] = i;
            }
        }

        // Unrolled level state transitions for edit distances 1..MaxDist
        #pragma unroll
        for (uint32_t distance = 1; distance <= MaxDist; ++distance) {
            previous_state[distance] = current_state[distance];

            // Bit-parallel DP operations: deletion, insertion, substitution, match
            uint64_t deletion = previous_state[distance - 1];
            uint64_t insertion = current_state[distance - 1] << 1;
            uint64_t substitution = previous_state[distance - 1] << 1;
            uint64_t match = (previous_state[distance] << 1) | current_pattern_mask;

            uint64_t new_state = deletion & insertion & substitution & match;

            if constexpr (Policy == MatchPolicy::kFirst) {
                // Latch on first occurrence at this edit distance
                if ((best_position[distance] == kCudaNotFound) && ((new_state & pattern.msb) == 0)) {
                    best_position[distance] = i;
                }
            } else {
                // Overwrite on every hit — keeps the latest (widest) position
                if ((new_state & pattern.msb) == 0) {
                    best_position[distance] = i;
                }
            }

            current_state[distance] = new_state;
        }
    }

    // Return best position at the lowest edit distance that matched
    constexpr uint32_t kStartDist = (Policy == MatchPolicy::kFirst) ? 1 : 0;
    #pragma unroll
    for (uint32_t distance = kStartDist; distance <= MaxDist; ++distance) {
        if (best_position[distance] != kCudaNotFound) {
            return (pattern.direction == SearchDirection::kBackward)
                ? (end - best_position[distance])
                : (begin + best_position[distance]);
        }
    }

    return kCudaNoMatchPosition;
}

/**
 * @brief Find both start and end positions of the best match in the reference.
 * Mirrors Bitap::FindStartEnd in bitap.cpp.
 */
template <uint32_t MaxDist>
__host__ __device__ __forceinline__ void bitap_find_start_end(
    const BitapPattern& fw_pattern,
    const BitapPattern& bw_pattern,
    const char* __restrict__ ref,
    int32_t begin,
    int32_t end,
    int32_t& out_start,
    int32_t& out_end
) {
    out_start = kCudaNoMatchPosition;
    out_end = kCudaNoMatchPosition;

    // Step 1: Forward pass locates the end boundary
    int32_t end_pos = bitap_search_window<MaxDist, MatchPolicy::kFirst>(fw_pattern, ref, begin, end);
    if (end_pos == kCudaNoMatchPosition) {
        return;
    }

    // Step 2: Backward pass recovers the start boundary within narrowed window
    const int32_t query_len = fw_pattern.query_len;
    const int32_t max_dist = static_cast<int32_t>(MaxDist);
    const int32_t search_begin = (end_pos - query_len - max_dist < begin) ? begin : (end_pos - query_len - max_dist);

    int32_t start_pos = bitap_search_window<MaxDist, MatchPolicy::kFirst>(bw_pattern, ref, search_begin, end_pos);
    if (start_pos == kCudaNoMatchPosition) {
        return;
    }

    out_start = start_pos;
    out_end = end_pos;
}

} // namespace xoos::demux::cuda
