// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/planning/proto/spiral_curve_config.proto

#include "modules/planning/proto/spiral_curve_config.pb.h"

#include <algorithm>

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/extension_set.h>
#include <google/protobuf/wire_format_lite.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>

PROTOBUF_PRAGMA_INIT_SEG

namespace _pb = ::PROTOBUF_NAMESPACE_ID;
namespace _pbi = _pb::internal;

namespace apollo {
namespace planning {
PROTOBUF_CONSTEXPR SpiralCurveConfig::SpiralCurveConfig(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_._has_bits_)*/{}
  , /*decltype(_impl_._cached_size_)*/{}
  , /*decltype(_impl_.newton_raphson_tol_)*/0.01
  , /*decltype(_impl_.simpson_size_)*/9
  , /*decltype(_impl_.newton_raphson_max_iter_)*/20} {}
struct SpiralCurveConfigDefaultTypeInternal {
  PROTOBUF_CONSTEXPR SpiralCurveConfigDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~SpiralCurveConfigDefaultTypeInternal() {}
  union {
    SpiralCurveConfig _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 SpiralCurveConfigDefaultTypeInternal _SpiralCurveConfig_default_instance_;
}  // namespace planning
}  // namespace apollo
static ::_pb::Metadata file_level_metadata_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto[1];
static constexpr ::_pb::EnumDescriptor const** file_level_enum_descriptors_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto = nullptr;
static constexpr ::_pb::ServiceDescriptor const** file_level_service_descriptors_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto = nullptr;

const uint32_t TableStruct_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::apollo::planning::SpiralCurveConfig, _impl_._has_bits_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::SpiralCurveConfig, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::apollo::planning::SpiralCurveConfig, _impl_.simpson_size_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::SpiralCurveConfig, _impl_.newton_raphson_tol_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::SpiralCurveConfig, _impl_.newton_raphson_max_iter_),
  1,
  0,
  2,
};
static const ::_pbi::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 9, -1, sizeof(::apollo::planning::SpiralCurveConfig)},
};

static const ::_pb::Message* const file_default_instances[] = {
  &::apollo::planning::_SpiralCurveConfig_default_instance_._instance,
};

const char descriptor_table_protodef_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n0modules/planning/proto/spiral_curve_co"
  "nfig.proto\022\017apollo.planning\"s\n\021SpiralCur"
  "veConfig\022\027\n\014simpson_size\030\001 \001(\005:\0019\022 \n\022new"
  "ton_raphson_tol\030\002 \001(\001:\0040.01\022#\n\027newton_ra"
  "phson_max_iter\030\003 \001(\005:\00220"
  ;
static ::_pbi::once_flag descriptor_table_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto_once;
const ::_pbi::DescriptorTable descriptor_table_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto = {
    false, false, 184, descriptor_table_protodef_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto,
    "modules/planning/proto/spiral_curve_config.proto",
    &descriptor_table_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto_once, nullptr, 0, 1,
    schemas, file_default_instances, TableStruct_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto::offsets,
    file_level_metadata_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto, file_level_enum_descriptors_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto,
    file_level_service_descriptors_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto,
};
PROTOBUF_ATTRIBUTE_WEAK const ::_pbi::DescriptorTable* descriptor_table_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto_getter() {
  return &descriptor_table_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto;
}

// Force running AddDescriptors() at dynamic initialization time.
PROTOBUF_ATTRIBUTE_INIT_PRIORITY2 static ::_pbi::AddDescriptorsRunner dynamic_init_dummy_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto(&descriptor_table_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto);
namespace apollo {
namespace planning {

// ===================================================================

class SpiralCurveConfig::_Internal {
 public:
  using HasBits = decltype(std::declval<SpiralCurveConfig>()._impl_._has_bits_);
  static void set_has_simpson_size(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static void set_has_newton_raphson_tol(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_newton_raphson_max_iter(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
};

SpiralCurveConfig::SpiralCurveConfig(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:apollo.planning.SpiralCurveConfig)
}
SpiralCurveConfig::SpiralCurveConfig(const SpiralCurveConfig& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){from._impl_._has_bits_}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.newton_raphson_tol_){}
    , decltype(_impl_.simpson_size_){}
    , decltype(_impl_.newton_raphson_max_iter_){}};

  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::memcpy(&_impl_.newton_raphson_tol_, &from._impl_.newton_raphson_tol_,
    static_cast<size_t>(reinterpret_cast<char*>(&_impl_.newton_raphson_max_iter_) -
    reinterpret_cast<char*>(&_impl_.newton_raphson_tol_)) + sizeof(_impl_.newton_raphson_max_iter_));
  // @@protoc_insertion_point(copy_constructor:apollo.planning.SpiralCurveConfig)
}

inline void SpiralCurveConfig::SharedCtor(
    ::_pb::Arena* arena, bool is_message_owned) {
  (void)arena;
  (void)is_message_owned;
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.newton_raphson_tol_){0.01}
    , decltype(_impl_.simpson_size_){9}
    , decltype(_impl_.newton_raphson_max_iter_){20}
  };
}

SpiralCurveConfig::~SpiralCurveConfig() {
  // @@protoc_insertion_point(destructor:apollo.planning.SpiralCurveConfig)
  if (auto *arena = _internal_metadata_.DeleteReturnArena<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>()) {
  (void)arena;
    return;
  }
  SharedDtor();
}

inline void SpiralCurveConfig::SharedDtor() {
  GOOGLE_DCHECK(GetArenaForAllocation() == nullptr);
}

void SpiralCurveConfig::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void SpiralCurveConfig::Clear() {
// @@protoc_insertion_point(message_clear_start:apollo.planning.SpiralCurveConfig)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x00000007u) {
    _impl_.newton_raphson_tol_ = 0.01;
    _impl_.simpson_size_ = 9;
    _impl_.newton_raphson_max_iter_ = 20;
  }
  _impl_._has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* SpiralCurveConfig::_InternalParse(const char* ptr, ::_pbi::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  while (!ctx->Done(&ptr)) {
    uint32_t tag;
    ptr = ::_pbi::ReadTag(ptr, &tag);
    switch (tag >> 3) {
      // optional int32 simpson_size = 1 [default = 9];
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 8)) {
          _Internal::set_has_simpson_size(&has_bits);
          _impl_.simpson_size_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint32(&ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional double newton_raphson_tol = 2 [default = 0.01];
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 17)) {
          _Internal::set_has_newton_raphson_tol(&has_bits);
          _impl_.newton_raphson_tol_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
        } else
          goto handle_unusual;
        continue;
      // optional int32 newton_raphson_max_iter = 3 [default = 20];
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 24)) {
          _Internal::set_has_newton_raphson_max_iter(&has_bits);
          _impl_.newton_raphson_max_iter_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint32(&ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      default:
        goto handle_unusual;
    }  // switch
  handle_unusual:
    if ((tag == 0) || ((tag & 7) == 4)) {
      CHK_(ptr);
      ctx->SetLastTag(tag);
      goto message_done;
    }
    ptr = UnknownFieldParse(
        tag,
        _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(),
        ptr, ctx);
    CHK_(ptr != nullptr);
  }  // while
message_done:
  _impl_._has_bits_.Or(has_bits);
  return ptr;
failure:
  ptr = nullptr;
  goto message_done;
#undef CHK_
}

uint8_t* SpiralCurveConfig::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:apollo.planning.SpiralCurveConfig)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  // optional int32 simpson_size = 1 [default = 9];
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteInt32ToArray(1, this->_internal_simpson_size(), target);
  }

  // optional double newton_raphson_tol = 2 [default = 0.01];
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteDoubleToArray(2, this->_internal_newton_raphson_tol(), target);
  }

  // optional int32 newton_raphson_max_iter = 3 [default = 20];
  if (cached_has_bits & 0x00000004u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteInt32ToArray(3, this->_internal_newton_raphson_max_iter(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:apollo.planning.SpiralCurveConfig)
  return target;
}

size_t SpiralCurveConfig::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:apollo.planning.SpiralCurveConfig)
  size_t total_size = 0;

  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x00000007u) {
    // optional double newton_raphson_tol = 2 [default = 0.01];
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 + 8;
    }

    // optional int32 simpson_size = 1 [default = 9];
    if (cached_has_bits & 0x00000002u) {
      total_size += ::_pbi::WireFormatLite::Int32SizePlusOne(this->_internal_simpson_size());
    }

    // optional int32 newton_raphson_max_iter = 3 [default = 20];
    if (cached_has_bits & 0x00000004u) {
      total_size += ::_pbi::WireFormatLite::Int32SizePlusOne(this->_internal_newton_raphson_max_iter());
    }

  }
  return MaybeComputeUnknownFieldsSize(total_size, &_impl_._cached_size_);
}

const ::PROTOBUF_NAMESPACE_ID::Message::ClassData SpiralCurveConfig::_class_data_ = {
    ::PROTOBUF_NAMESPACE_ID::Message::CopyWithSizeCheck,
    SpiralCurveConfig::MergeImpl
};
const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*SpiralCurveConfig::GetClassData() const { return &_class_data_; }

void SpiralCurveConfig::MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to,
                      const ::PROTOBUF_NAMESPACE_ID::Message& from) {
  static_cast<SpiralCurveConfig *>(to)->MergeFrom(
      static_cast<const SpiralCurveConfig &>(from));
}


void SpiralCurveConfig::MergeFrom(const SpiralCurveConfig& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:apollo.planning.SpiralCurveConfig)
  GOOGLE_DCHECK_NE(&from, this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._impl_._has_bits_[0];
  if (cached_has_bits & 0x00000007u) {
    if (cached_has_bits & 0x00000001u) {
      _impl_.newton_raphson_tol_ = from._impl_.newton_raphson_tol_;
    }
    if (cached_has_bits & 0x00000002u) {
      _impl_.simpson_size_ = from._impl_.simpson_size_;
    }
    if (cached_has_bits & 0x00000004u) {
      _impl_.newton_raphson_max_iter_ = from._impl_.newton_raphson_max_iter_;
    }
    _impl_._has_bits_[0] |= cached_has_bits;
  }
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void SpiralCurveConfig::CopyFrom(const SpiralCurveConfig& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:apollo.planning.SpiralCurveConfig)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool SpiralCurveConfig::IsInitialized() const {
  return true;
}

void SpiralCurveConfig::InternalSwap(SpiralCurveConfig* other) {
  using std::swap;
  _internal_metadata_.InternalSwap(&other->_internal_metadata_);
  swap(_impl_._has_bits_[0], other->_impl_._has_bits_[0]);
  swap(_impl_.newton_raphson_tol_, other->_impl_.newton_raphson_tol_);
  swap(_impl_.simpson_size_, other->_impl_.simpson_size_);
  swap(_impl_.newton_raphson_max_iter_, other->_impl_.newton_raphson_max_iter_);
}

::PROTOBUF_NAMESPACE_ID::Metadata SpiralCurveConfig::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto_getter, &descriptor_table_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto_once,
      file_level_metadata_modules_2fplanning_2fproto_2fspiral_5fcurve_5fconfig_2eproto[0]);
}

// @@protoc_insertion_point(namespace_scope)
}  // namespace planning
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::apollo::planning::SpiralCurveConfig*
Arena::CreateMaybeMessage< ::apollo::planning::SpiralCurveConfig >(Arena* arena) {
  return Arena::CreateMessageInternal< ::apollo::planning::SpiralCurveConfig >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
