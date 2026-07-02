#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json_fwd.hpp>

#include "xoos/types/int.h"

namespace fs = std::filesystem;

namespace xoos::demux {
enum class BarcodeType {
  kRunway,
  kUmi,
  kSid,
  kStem,
  kAnchor,
  kLoop,
};

std::string Format(BarcodeType type);

enum class LutTransform {
  kReverse,
  kReverseComplement,
  kNone,
};

/**
 * Define how the LUT is generated from a set of input sequences.
 *
 * Exactly one of `sequences` (file path) or `sequence` (inline string) must be set.
 * Use `sequence` for fixed single-record adapter components (runway, stem, anchor, loop).
 * Use `sequences` for multi-record barcode pools (SID, UMI).
 */
struct LutDefinition {
  /// @brief A FASTA file containing all of the sequences and their names for this LUT
  fs::path sequences;
  /// @brief A single inline DNA sequence, used instead of a FASTA file for fixed adapter components
  std::string sequence;
  /// @brief The LUT will contain all sequences within this max_edit_distance of the sequences
  u32 max_edit_distance{0};
  /// @brief The sequences can first be transformed before generating the LUT
  LutTransform transform{LutTransform::kNone};
  /// @brief Bases prepended to each FASTA sequence before transform and LUT generation
  std::string prefix;
  /// @brief Bases appended to each FASTA sequence before transform and LUT generation
  std::string suffix;

  /// @brief True when the LUT source is an inline sequence rather than a FASTA file
  bool HasInlineSequence() const { return !sequence.empty(); }

  friend bool operator==(const LutDefinition& a, const LutDefinition& b) = default;
};

/**
 * Define how to search for a barcode in an adapter, search to the left and right of the expected position by the
 * defined amount.
 */
struct SearchDefinition {
  int max_wiggle_left;
  int max_wiggle_right;

  friend bool operator==(const SearchDefinition& a, const SearchDefinition& b) = default;
};

struct BarcodeDefinition {
  BarcodeType type{};
  LutDefinition lut;
  SearchDefinition search{};

  friend bool operator==(const BarcodeDefinition& a, const BarcodeDefinition& b) = default;
};

std::string FormatJson(const BarcodeDefinition& bd);

enum class AdapterType { kDuplex, kDuplexUMI, kDuplexStem, kYsu, kYs, kSimplex, kSimplex10x };

struct AdapterDesign {
  std::string name;
  AdapterType type{};
  std::vector<BarcodeDefinition> adapter_5p;
  std::vector<BarcodeDefinition> adapter_3p;
  std::optional<std::string> bait_5p;
  std::optional<std::string> bait_3p;
};

/**
 * A struct used to describe multiple potential adapter designs.
 * This struct will most often be loaded from a bundle containing multiple adapter designs,
 * it will describe the details of the demultiplexing and trimming logic used for that adapter design.
 */
struct AdapterDesignManifest {
  std::string default_adapter_design_name;
  std::vector<AdapterDesign> adapter_designs;
};

/**
 * Load adapter design with the given name from the given bundle, or use default if name is nullopt.
 *
 * Bundle must contain a "manifest.json" which contains an AdapterDesignManifest in JSON representation.
 */
AdapterDesign LoadAdapterDesign(const fs::path& adapter_design_bundle,
                                const std::optional<std::string>& adapter_design_name);

// NOLINTBEGIN(readability-identifier-naming) -- snake_case names required by nlohmann ADL contract
void to_json(nlohmann::json& j, const BarcodeType& e);          // NOSONAR(S100)
void from_json(const nlohmann::json& j, BarcodeType& e);        // NOSONAR(S100)
void to_json(nlohmann::json& j, const LutTransform& e);         // NOSONAR(S100)
void from_json(const nlohmann::json& j, LutTransform& e);       // NOSONAR(S100)
void to_json(nlohmann::json& j, const AdapterType& e);          // NOSONAR(S100)
void from_json(const nlohmann::json& j, AdapterType& e);        // NOSONAR(S100)
void to_json(nlohmann::json& j, const LutDefinition& d);        // NOSONAR(S100)
void from_json(const nlohmann::json& j, LutDefinition& d);      // NOSONAR(S100)
void to_json(nlohmann::json& j, const SearchDefinition& d);     // NOSONAR(S100)
void from_json(const nlohmann::json& j, SearchDefinition& d);   // NOSONAR(S100)
void to_json(nlohmann::json& j, const BarcodeDefinition& d);    // NOSONAR(S100)
void from_json(const nlohmann::json& j, BarcodeDefinition& d);  // NOSONAR(S100)
void to_json(nlohmann::json& j, const AdapterDesign& d);        // NOSONAR(S100)
void from_json(const nlohmann::json& j, AdapterDesign& d);      // NOSONAR(S100)
// NOLINTEND(readability-identifier-naming)

}  // namespace xoos::demux
