#include "adapter-design/adapter-design.h"

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/vfs/vfs-iterator.h>
#include <xoos/vfs/vfs.h>

#include <array>
#include <filesystem>
#include <format>
#include <memory>
#include <unordered_map>

#include <nlohmann/json.hpp>

#include "adapter-design/adapter-design-constants.h"

namespace vfs = xoos::vfs;

using json = nlohmann::json;  // NOLINT

namespace {

/// Length of the trailing separator ", " appended after each enum value name.
constexpr std::size_t kSeparatorLength = 2;

/// Lookup a JSON string value in a mapping table and return the corresponding enum.
/// Throws if the value is not found, listing all valid options.
template <typename EnumType, std::size_t N>
EnumType StrictEnumLookup(const json& j, const std::array<std::pair<EnumType, std::string_view>, N>& mapping,
                          std::string_view type_name) {
  const auto str = j.get<std::string>();
  for (const auto& [value, name] : mapping) {
    if (name == str) {
      return value;
    }
  }
  std::string valid;
  for (const auto& [_, name] : mapping) {
    valid += std::format("'{}', ", name);
  }
  if (valid.size() >= kSeparatorLength) {
    valid.resize(valid.size() - kSeparatorLength);
  }
  throw xoos::error::Error("Unknown {} '{}', valid values: {}", type_name, str, valid);
}

/// Lookup an enum value in a mapping table and return the corresponding JSON string.
/// Throws if the value is not found.
template <typename EnumType, std::size_t N>
std::string_view StrictEnumToString(EnumType e, const std::array<std::pair<EnumType, std::string_view>, N>& mapping,
                                    std::string_view type_name) {
  for (const auto& [value, name] : mapping) {
    if (value == e) {
      return name;
    }
  }
  throw xoos::error::Error("Unknown {} enum value {}", type_name, static_cast<xoos::s32>(e));
}

}  // namespace

namespace xoos::demux {

// --- Enum mapping tables ---

static constexpr auto kBarcodeTypeMapping = std::to_array<std::pair<BarcodeType, std::string_view>>({
    {BarcodeType::kLoop, "loop"},
    {BarcodeType::kRunway, "runway"},
    {BarcodeType::kUmi, "umi"},
    {BarcodeType::kSid, "sid"},
    {BarcodeType::kStem, "stem"},
    {BarcodeType::kAnchor, "anchor"},
});

static constexpr auto kLutTransformMapping = std::to_array<std::pair<LutTransform, std::string_view>>({
    {LutTransform::kNone, "none"},
    {LutTransform::kReverse, "reverse"},
    {LutTransform::kReverseComplement, "reverse_complement"},
});

static constexpr auto kAdapterTypeMapping = std::to_array<std::pair<AdapterType, std::string_view>>({
    {AdapterType::kDuplex, "duplex"},
    {AdapterType::kDuplexUMI, "duplex_umi"},
    {AdapterType::kDuplexStem, "duplex_stem"},
    {AdapterType::kYsu, "YSU"},
    {AdapterType::kYs, "YS"},
    {AdapterType::kSimplex, "simplex"},
    {AdapterType::kSimplex10x, "simplex-10x"},
});

// --- Enum to_json / from_json ---
// NOLINTBEGIN(readability-identifier-naming) -- snake_case names required by nlohmann ADL contract

void to_json(json& j, const BarcodeType& e) {  // NOSONAR(S100)
  j = StrictEnumToString(e, kBarcodeTypeMapping, "BarcodeType");
}

void from_json(const json& j, BarcodeType& e) {  // NOSONAR(S100)
  e = StrictEnumLookup(j, kBarcodeTypeMapping, "BarcodeType");
}

void to_json(json& j, const LutTransform& e) {  // NOSONAR(S100)
  j = StrictEnumToString(e, kLutTransformMapping, "LutTransform");
}

void from_json(const json& j, LutTransform& e) {  // NOSONAR(S100)
  e = StrictEnumLookup(j, kLutTransformMapping, "LutTransform");
}

void to_json(json& j, const AdapterType& e) {  // NOSONAR(S100)
  j = StrictEnumToString(e, kAdapterTypeMapping, "AdapterType");
}

void from_json(const json& j, AdapterType& e) {  // NOSONAR(S100)
  e = StrictEnumLookup(j, kAdapterTypeMapping, "AdapterType");
}

// --- Struct to_json / from_json ---

/**
 * @brief Custom JSON serialization for LutDefinition.
 *
 * Writes either "sequence" (inline string) or "sequences" (file path), plus optional
 * "prefix" and "suffix" fields.
 *
 * @param[out] j   JSON object to populate.
 * @param      d   LutDefinition to serialize.
 */
void to_json(json& j, const LutDefinition& d) {  // NOSONAR(S100)
  j = json{{kMaxEditDistanceKey, d.max_edit_distance}, {kTransformKey, d.transform}};
  if (d.HasInlineSequence()) {
    j[kSequenceKey] = d.sequence;
  } else {
    j[kSequencesKey] = d.sequences;
  }
  if (!d.prefix.empty()) {
    j[kPrefixKey] = d.prefix;
  }
  if (!d.suffix.empty()) {
    j[kSuffixKey] = d.suffix;
  }
}

/**
 * @brief Custom JSON deserialization for LutDefinition.
 *
 * Reads exactly one of "sequence" (inline string) or "sequences" (file path).
 * Errors if both or neither are present.
 *
 * @param      j   JSON object to read from.
 * @param[out] d   LutDefinition to populate.
 */
void from_json(const json& j, LutDefinition& d) {  // NOSONAR(S100)
  const bool has_sequence = j.contains(kSequenceKey);
  const bool has_sequences = j.contains(kSequencesKey);
  if (has_sequence == has_sequences) {
    throw error::Error("LutDefinition must have exactly one of '{}' or '{}', found {}", kSequenceKey, kSequencesKey,
                       has_sequence ? "both" : "neither");
  }
  if (has_sequence) {
    j.at(kSequenceKey).get_to(d.sequence);
    if (d.sequence.empty()) {
      throw error::Error("LutDefinition '{}' must not be empty", kSequenceKey);
    }
    d.sequences.clear();
  } else {
    j.at(kSequencesKey).get_to(d.sequences);
    if (d.sequences.empty()) {
      throw error::Error("LutDefinition '{}' must not be empty", kSequencesKey);
    }
    d.sequence.clear();
  }
  j.at(kMaxEditDistanceKey).get_to(d.max_edit_distance);
  j.at(kTransformKey).get_to(d.transform);
  d.prefix = j.value(kPrefixKey, std::string{});
  d.suffix = j.value(kSuffixKey, std::string{});
}

void to_json(json& j, const SearchDefinition& d) {  // NOSONAR(S100)
  j = json{{"max_wiggle_left", d.max_wiggle_left}, {"max_wiggle_right", d.max_wiggle_right}};
}

void from_json(const json& j, SearchDefinition& d) {  // NOSONAR(S100)
  j.at("max_wiggle_left").get_to(d.max_wiggle_left);
  j.at("max_wiggle_right").get_to(d.max_wiggle_right);
}

void to_json(json& j, const BarcodeDefinition& d) {  // NOSONAR(S100)
  j = json{{"type", d.type}, {"lut", d.lut}, {"search", d.search}};
}

void from_json(const json& j, BarcodeDefinition& d) {  // NOSONAR(S100)
  j.at("type").get_to(d.type);
  j.at("lut").get_to(d.lut);
  j.at("search").get_to(d.search);
}

void to_json(json& j, const AdapterDesign& d) {  // NOSONAR(S100)
  j = json{{kNameKey, d.name}, {kTypeKey, d.type}, {kAdapter5pKey, d.adapter_5p}, {kAdapter3pKey, d.adapter_3p}};
  if (d.bait_5p.has_value()) {
    j[kBait5pKey] = d.bait_5p.value();
  }
  if (d.bait_3p.has_value()) {
    j[kBait3pKey] = d.bait_3p.value();
  }
}

void from_json(const json& j, AdapterDesign& d) {  // NOSONAR(S100)
  j.at(kNameKey).get_to(d.name);
  j.at(kTypeKey).get_to(d.type);
  j.at(kAdapter5pKey).get_to(d.adapter_5p);
  j.at(kAdapter3pKey).get_to(d.adapter_3p);
  if (j.contains(kBait5pKey)) {
    d.bait_5p = j.at(kBait5pKey).get<std::string>();
  }
  if (j.contains(kBait3pKey)) {
    d.bait_3p = j.at(kBait3pKey).get<std::string>();
  }
}

void to_json(json& j, const AdapterDesignManifest& d) {  // NOSONAR(S100)
  j = json{{"default_adapter_design_name", d.default_adapter_design_name}, {"adapter_designs", d.adapter_designs}};
}

void from_json(const json& j, AdapterDesignManifest& d) {  // NOSONAR(S100)
  j.at("default_adapter_design_name").get_to(d.default_adapter_design_name);
  j.at("adapter_designs").get_to(d.adapter_designs);
}

// NOLINTEND(readability-identifier-naming)

std::string Format(BarcodeType type) {
  std::unordered_map<BarcodeType, std::string> barcode_type_names = {
      {BarcodeType::kRunway, "Runway"}, {BarcodeType::kSid, "Sid"},       {BarcodeType::kUmi, "Umi"},
      {BarcodeType::kStem, "Stem"},     {BarcodeType::kAnchor, "Anchor"}, {BarcodeType::kLoop, "Loop"}};
  return barcode_type_names.at(type);
}

std::string FormatJson(const BarcodeDefinition& bd) {
  json j;
  to_json(j, bd);
  return j.dump();
}

AdapterDesignManifest LoadBundleManifest(const vfs::VirtualFilesystemPtr& vfs, const fs::path& path);

AdapterDesign FindAdapterDesign(const AdapterDesignManifest& manifest, const std::string& design_name);

AdapterDesign LoadAdapterDesign(const fs::path& adapter_design_bundle,
                                const std::optional<std::string>& adapter_design_name) {
  vfs::VirtualFilesystemPtr vfs = vfs::Open(adapter_design_bundle);
  AdapterDesignManifest manifest = LoadBundleManifest(vfs, fs::path{"manifest.json"});
  std::string design_name = adapter_design_name.value_or(manifest.default_adapter_design_name);
  Logging::Info("Loading adapter design '{}' from bundle '{}'", design_name, adapter_design_bundle.string());
  return FindAdapterDesign(manifest, design_name);
}

/**
 * Load JSON manifest from bundle.
 */
AdapterDesignManifest LoadBundleManifest(const vfs::VirtualFilesystemPtr& vfs, const fs::path& path) {
  vfs::VirtualFileHandlePtr fh = vfs->Open(path);
  if (fh == nullptr) {
    throw error::Error("Manifest '{}' does not exist in bundle '{}'", path.string(), vfs->GetName());
  }

  AdapterDesignManifest manifest;
  try {
    from_json(json::parse(std::make_shared<vfs::VirtualFileStream>(fh)), manifest);
  } catch (const json::exception& ex) {
    throw error::Error("Manifest '{}' in bundle '{}' cannot be parsed due to '{}'", path.string(), vfs->GetName(),
                       ex.what());
  }
  return manifest;
}

/**
 * Find adapter design with provided name, if there are multiple adapter design with the same name
 * only the first will be found.
 */
AdapterDesign FindAdapterDesign(const AdapterDesignManifest& manifest, const std::string& design_name) {
  auto adapter_design = std::ranges::find_if(manifest.adapter_designs,
                                             [&design_name](const auto& item) { return item.name == design_name; });
  if (adapter_design == manifest.adapter_designs.end()) {
    std::string available_designs;
    for (const auto& design : manifest.adapter_designs) {
      available_designs += std::format(" {} type: {}", design.name, static_cast<u32>(design.type));
    }
    throw error::Error("Adapter design '{}' does not exist in manifest, available designs:{}", design_name,
                       available_designs);
  }
  return *adapter_design;
}
}  // namespace xoos::demux
