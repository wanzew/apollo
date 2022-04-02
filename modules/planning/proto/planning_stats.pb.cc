// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/planning/proto/planning_stats.proto

#include "modules/planning/proto/planning_stats.pb.h"

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
PROTOBUF_CONSTEXPR StatsGroup::StatsGroup(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_._has_bits_)*/{}
  , /*decltype(_impl_._cached_size_)*/{}
  , /*decltype(_impl_.max_)*/0
  , /*decltype(_impl_.sum_)*/0
  , /*decltype(_impl_.avg_)*/0
  , /*decltype(_impl_.num_)*/0
  , /*decltype(_impl_.min_)*/10000000000} {}
struct StatsGroupDefaultTypeInternal {
  PROTOBUF_CONSTEXPR StatsGroupDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~StatsGroupDefaultTypeInternal() {}
  union {
    StatsGroup _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 StatsGroupDefaultTypeInternal _StatsGroup_default_instance_;
PROTOBUF_CONSTEXPR PlanningStats::PlanningStats(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_._has_bits_)*/{}
  , /*decltype(_impl_._cached_size_)*/{}
  , /*decltype(_impl_.total_path_length_)*/nullptr
  , /*decltype(_impl_.total_path_time_)*/nullptr
  , /*decltype(_impl_.v_)*/nullptr
  , /*decltype(_impl_.a_)*/nullptr
  , /*decltype(_impl_.kappa_)*/nullptr
  , /*decltype(_impl_.dkappa_)*/nullptr} {}
struct PlanningStatsDefaultTypeInternal {
  PROTOBUF_CONSTEXPR PlanningStatsDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~PlanningStatsDefaultTypeInternal() {}
  union {
    PlanningStats _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 PlanningStatsDefaultTypeInternal _PlanningStats_default_instance_;
}  // namespace planning
}  // namespace apollo
static ::_pb::Metadata file_level_metadata_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto[2];
static constexpr ::_pb::EnumDescriptor const** file_level_enum_descriptors_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto = nullptr;
static constexpr ::_pb::ServiceDescriptor const** file_level_service_descriptors_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto = nullptr;

const uint32_t TableStruct_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::apollo::planning::StatsGroup, _impl_._has_bits_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::StatsGroup, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::apollo::planning::StatsGroup, _impl_.max_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::StatsGroup, _impl_.min_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::StatsGroup, _impl_.sum_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::StatsGroup, _impl_.avg_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::StatsGroup, _impl_.num_),
  0,
  4,
  1,
  2,
  3,
  PROTOBUF_FIELD_OFFSET(::apollo::planning::PlanningStats, _impl_._has_bits_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::PlanningStats, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::apollo::planning::PlanningStats, _impl_.total_path_length_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::PlanningStats, _impl_.total_path_time_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::PlanningStats, _impl_.v_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::PlanningStats, _impl_.a_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::PlanningStats, _impl_.kappa_),
  PROTOBUF_FIELD_OFFSET(::apollo::planning::PlanningStats, _impl_.dkappa_),
  0,
  1,
  2,
  3,
  4,
  5,
};
static const ::_pbi::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 11, -1, sizeof(::apollo::planning::StatsGroup)},
  { 16, 28, -1, sizeof(::apollo::planning::PlanningStats)},
};

static const ::_pb::Message* const file_default_instances[] = {
  &::apollo::planning::_StatsGroup_default_instance_._instance,
  &::apollo::planning::_PlanningStats_default_instance_._instance,
};

const char descriptor_table_protodef_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n+modules/planning/proto/planning_stats."
  "proto\022\017apollo.planning\"Z\n\nStatsGroup\022\013\n\003"
  "max\030\001 \001(\001\022\030\n\003min\030\002 \001(\001:\01310000000000\022\013\n\003s"
  "um\030\003 \001(\001\022\013\n\003avg\030\004 \001(\001\022\013\n\003num\030\005 \001(\005\"\246\002\n\rP"
  "lanningStats\0226\n\021total_path_length\030\001 \001(\0132"
  "\033.apollo.planning.StatsGroup\0224\n\017total_pa"
  "th_time\030\002 \001(\0132\033.apollo.planning.StatsGro"
  "up\022&\n\001v\030\003 \001(\0132\033.apollo.planning.StatsGro"
  "up\022&\n\001a\030\004 \001(\0132\033.apollo.planning.StatsGro"
  "up\022*\n\005kappa\030\005 \001(\0132\033.apollo.planning.Stat"
  "sGroup\022+\n\006dkappa\030\006 \001(\0132\033.apollo.planning"
  ".StatsGroup"
  ;
static ::_pbi::once_flag descriptor_table_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto_once;
const ::_pbi::DescriptorTable descriptor_table_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto = {
    false, false, 451, descriptor_table_protodef_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto,
    "modules/planning/proto/planning_stats.proto",
    &descriptor_table_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto_once, nullptr, 0, 2,
    schemas, file_default_instances, TableStruct_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto::offsets,
    file_level_metadata_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto, file_level_enum_descriptors_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto,
    file_level_service_descriptors_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto,
};
PROTOBUF_ATTRIBUTE_WEAK const ::_pbi::DescriptorTable* descriptor_table_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto_getter() {
  return &descriptor_table_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto;
}

// Force running AddDescriptors() at dynamic initialization time.
PROTOBUF_ATTRIBUTE_INIT_PRIORITY2 static ::_pbi::AddDescriptorsRunner dynamic_init_dummy_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto(&descriptor_table_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto);
namespace apollo {
namespace planning {

// ===================================================================

class StatsGroup::_Internal {
 public:
  using HasBits = decltype(std::declval<StatsGroup>()._impl_._has_bits_);
  static void set_has_max(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_min(HasBits* has_bits) {
    (*has_bits)[0] |= 16u;
  }
  static void set_has_sum(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static void set_has_avg(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static void set_has_num(HasBits* has_bits) {
    (*has_bits)[0] |= 8u;
  }
};

StatsGroup::StatsGroup(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:apollo.planning.StatsGroup)
}
StatsGroup::StatsGroup(const StatsGroup& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){from._impl_._has_bits_}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.max_){}
    , decltype(_impl_.sum_){}
    , decltype(_impl_.avg_){}
    , decltype(_impl_.num_){}
    , decltype(_impl_.min_){}};

  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::memcpy(&_impl_.max_, &from._impl_.max_,
    static_cast<size_t>(reinterpret_cast<char*>(&_impl_.min_) -
    reinterpret_cast<char*>(&_impl_.max_)) + sizeof(_impl_.min_));
  // @@protoc_insertion_point(copy_constructor:apollo.planning.StatsGroup)
}

inline void StatsGroup::SharedCtor(
    ::_pb::Arena* arena, bool is_message_owned) {
  (void)arena;
  (void)is_message_owned;
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.max_){0}
    , decltype(_impl_.sum_){0}
    , decltype(_impl_.avg_){0}
    , decltype(_impl_.num_){0}
    , decltype(_impl_.min_){10000000000}
  };
}

StatsGroup::~StatsGroup() {
  // @@protoc_insertion_point(destructor:apollo.planning.StatsGroup)
  if (auto *arena = _internal_metadata_.DeleteReturnArena<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>()) {
  (void)arena;
    return;
  }
  SharedDtor();
}

inline void StatsGroup::SharedDtor() {
  GOOGLE_DCHECK(GetArenaForAllocation() == nullptr);
}

void StatsGroup::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void StatsGroup::Clear() {
// @@protoc_insertion_point(message_clear_start:apollo.planning.StatsGroup)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    ::memset(&_impl_.max_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&_impl_.num_) -
        reinterpret_cast<char*>(&_impl_.max_)) + sizeof(_impl_.num_));
    _impl_.min_ = 10000000000;
  }
  _impl_._has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* StatsGroup::_InternalParse(const char* ptr, ::_pbi::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  while (!ctx->Done(&ptr)) {
    uint32_t tag;
    ptr = ::_pbi::ReadTag(ptr, &tag);
    switch (tag >> 3) {
      // optional double max = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 9)) {
          _Internal::set_has_max(&has_bits);
          _impl_.max_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
        } else
          goto handle_unusual;
        continue;
      // optional double min = 2 [default = 10000000000];
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 17)) {
          _Internal::set_has_min(&has_bits);
          _impl_.min_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
        } else
          goto handle_unusual;
        continue;
      // optional double sum = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 25)) {
          _Internal::set_has_sum(&has_bits);
          _impl_.sum_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
        } else
          goto handle_unusual;
        continue;
      // optional double avg = 4;
      case 4:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 33)) {
          _Internal::set_has_avg(&has_bits);
          _impl_.avg_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
        } else
          goto handle_unusual;
        continue;
      // optional int32 num = 5;
      case 5:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 40)) {
          _Internal::set_has_num(&has_bits);
          _impl_.num_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint32(&ptr);
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

uint8_t* StatsGroup::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:apollo.planning.StatsGroup)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  // optional double max = 1;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteDoubleToArray(1, this->_internal_max(), target);
  }

  // optional double min = 2 [default = 10000000000];
  if (cached_has_bits & 0x00000010u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteDoubleToArray(2, this->_internal_min(), target);
  }

  // optional double sum = 3;
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteDoubleToArray(3, this->_internal_sum(), target);
  }

  // optional double avg = 4;
  if (cached_has_bits & 0x00000004u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteDoubleToArray(4, this->_internal_avg(), target);
  }

  // optional int32 num = 5;
  if (cached_has_bits & 0x00000008u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteInt32ToArray(5, this->_internal_num(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:apollo.planning.StatsGroup)
  return target;
}

size_t StatsGroup::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:apollo.planning.StatsGroup)
  size_t total_size = 0;

  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    // optional double max = 1;
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 + 8;
    }

    // optional double sum = 3;
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 + 8;
    }

    // optional double avg = 4;
    if (cached_has_bits & 0x00000004u) {
      total_size += 1 + 8;
    }

    // optional int32 num = 5;
    if (cached_has_bits & 0x00000008u) {
      total_size += ::_pbi::WireFormatLite::Int32SizePlusOne(this->_internal_num());
    }

    // optional double min = 2 [default = 10000000000];
    if (cached_has_bits & 0x00000010u) {
      total_size += 1 + 8;
    }

  }
  return MaybeComputeUnknownFieldsSize(total_size, &_impl_._cached_size_);
}

const ::PROTOBUF_NAMESPACE_ID::Message::ClassData StatsGroup::_class_data_ = {
    ::PROTOBUF_NAMESPACE_ID::Message::CopyWithSizeCheck,
    StatsGroup::MergeImpl
};
const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*StatsGroup::GetClassData() const { return &_class_data_; }

void StatsGroup::MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to,
                      const ::PROTOBUF_NAMESPACE_ID::Message& from) {
  static_cast<StatsGroup *>(to)->MergeFrom(
      static_cast<const StatsGroup &>(from));
}


void StatsGroup::MergeFrom(const StatsGroup& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:apollo.planning.StatsGroup)
  GOOGLE_DCHECK_NE(&from, this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._impl_._has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    if (cached_has_bits & 0x00000001u) {
      _impl_.max_ = from._impl_.max_;
    }
    if (cached_has_bits & 0x00000002u) {
      _impl_.sum_ = from._impl_.sum_;
    }
    if (cached_has_bits & 0x00000004u) {
      _impl_.avg_ = from._impl_.avg_;
    }
    if (cached_has_bits & 0x00000008u) {
      _impl_.num_ = from._impl_.num_;
    }
    if (cached_has_bits & 0x00000010u) {
      _impl_.min_ = from._impl_.min_;
    }
    _impl_._has_bits_[0] |= cached_has_bits;
  }
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void StatsGroup::CopyFrom(const StatsGroup& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:apollo.planning.StatsGroup)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool StatsGroup::IsInitialized() const {
  return true;
}

void StatsGroup::InternalSwap(StatsGroup* other) {
  using std::swap;
  _internal_metadata_.InternalSwap(&other->_internal_metadata_);
  swap(_impl_._has_bits_[0], other->_impl_._has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(StatsGroup, _impl_.num_)
      + sizeof(StatsGroup::_impl_.num_)
      - PROTOBUF_FIELD_OFFSET(StatsGroup, _impl_.max_)>(
          reinterpret_cast<char*>(&_impl_.max_),
          reinterpret_cast<char*>(&other->_impl_.max_));
  swap(_impl_.min_, other->_impl_.min_);
}

::PROTOBUF_NAMESPACE_ID::Metadata StatsGroup::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto_getter, &descriptor_table_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto_once,
      file_level_metadata_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto[0]);
}

// ===================================================================

class PlanningStats::_Internal {
 public:
  using HasBits = decltype(std::declval<PlanningStats>()._impl_._has_bits_);
  static const ::apollo::planning::StatsGroup& total_path_length(const PlanningStats* msg);
  static void set_has_total_path_length(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static const ::apollo::planning::StatsGroup& total_path_time(const PlanningStats* msg);
  static void set_has_total_path_time(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static const ::apollo::planning::StatsGroup& v(const PlanningStats* msg);
  static void set_has_v(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static const ::apollo::planning::StatsGroup& a(const PlanningStats* msg);
  static void set_has_a(HasBits* has_bits) {
    (*has_bits)[0] |= 8u;
  }
  static const ::apollo::planning::StatsGroup& kappa(const PlanningStats* msg);
  static void set_has_kappa(HasBits* has_bits) {
    (*has_bits)[0] |= 16u;
  }
  static const ::apollo::planning::StatsGroup& dkappa(const PlanningStats* msg);
  static void set_has_dkappa(HasBits* has_bits) {
    (*has_bits)[0] |= 32u;
  }
};

const ::apollo::planning::StatsGroup&
PlanningStats::_Internal::total_path_length(const PlanningStats* msg) {
  return *msg->_impl_.total_path_length_;
}
const ::apollo::planning::StatsGroup&
PlanningStats::_Internal::total_path_time(const PlanningStats* msg) {
  return *msg->_impl_.total_path_time_;
}
const ::apollo::planning::StatsGroup&
PlanningStats::_Internal::v(const PlanningStats* msg) {
  return *msg->_impl_.v_;
}
const ::apollo::planning::StatsGroup&
PlanningStats::_Internal::a(const PlanningStats* msg) {
  return *msg->_impl_.a_;
}
const ::apollo::planning::StatsGroup&
PlanningStats::_Internal::kappa(const PlanningStats* msg) {
  return *msg->_impl_.kappa_;
}
const ::apollo::planning::StatsGroup&
PlanningStats::_Internal::dkappa(const PlanningStats* msg) {
  return *msg->_impl_.dkappa_;
}
PlanningStats::PlanningStats(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:apollo.planning.PlanningStats)
}
PlanningStats::PlanningStats(const PlanningStats& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){from._impl_._has_bits_}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.total_path_length_){nullptr}
    , decltype(_impl_.total_path_time_){nullptr}
    , decltype(_impl_.v_){nullptr}
    , decltype(_impl_.a_){nullptr}
    , decltype(_impl_.kappa_){nullptr}
    , decltype(_impl_.dkappa_){nullptr}};

  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  if (from._internal_has_total_path_length()) {
    _impl_.total_path_length_ = new ::apollo::planning::StatsGroup(*from._impl_.total_path_length_);
  }
  if (from._internal_has_total_path_time()) {
    _impl_.total_path_time_ = new ::apollo::planning::StatsGroup(*from._impl_.total_path_time_);
  }
  if (from._internal_has_v()) {
    _impl_.v_ = new ::apollo::planning::StatsGroup(*from._impl_.v_);
  }
  if (from._internal_has_a()) {
    _impl_.a_ = new ::apollo::planning::StatsGroup(*from._impl_.a_);
  }
  if (from._internal_has_kappa()) {
    _impl_.kappa_ = new ::apollo::planning::StatsGroup(*from._impl_.kappa_);
  }
  if (from._internal_has_dkappa()) {
    _impl_.dkappa_ = new ::apollo::planning::StatsGroup(*from._impl_.dkappa_);
  }
  // @@protoc_insertion_point(copy_constructor:apollo.planning.PlanningStats)
}

inline void PlanningStats::SharedCtor(
    ::_pb::Arena* arena, bool is_message_owned) {
  (void)arena;
  (void)is_message_owned;
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.total_path_length_){nullptr}
    , decltype(_impl_.total_path_time_){nullptr}
    , decltype(_impl_.v_){nullptr}
    , decltype(_impl_.a_){nullptr}
    , decltype(_impl_.kappa_){nullptr}
    , decltype(_impl_.dkappa_){nullptr}
  };
}

PlanningStats::~PlanningStats() {
  // @@protoc_insertion_point(destructor:apollo.planning.PlanningStats)
  if (auto *arena = _internal_metadata_.DeleteReturnArena<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>()) {
  (void)arena;
    return;
  }
  SharedDtor();
}

inline void PlanningStats::SharedDtor() {
  GOOGLE_DCHECK(GetArenaForAllocation() == nullptr);
  if (this != internal_default_instance()) delete _impl_.total_path_length_;
  if (this != internal_default_instance()) delete _impl_.total_path_time_;
  if (this != internal_default_instance()) delete _impl_.v_;
  if (this != internal_default_instance()) delete _impl_.a_;
  if (this != internal_default_instance()) delete _impl_.kappa_;
  if (this != internal_default_instance()) delete _impl_.dkappa_;
}

void PlanningStats::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void PlanningStats::Clear() {
// @@protoc_insertion_point(message_clear_start:apollo.planning.PlanningStats)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x0000003fu) {
    if (cached_has_bits & 0x00000001u) {
      GOOGLE_DCHECK(_impl_.total_path_length_ != nullptr);
      _impl_.total_path_length_->Clear();
    }
    if (cached_has_bits & 0x00000002u) {
      GOOGLE_DCHECK(_impl_.total_path_time_ != nullptr);
      _impl_.total_path_time_->Clear();
    }
    if (cached_has_bits & 0x00000004u) {
      GOOGLE_DCHECK(_impl_.v_ != nullptr);
      _impl_.v_->Clear();
    }
    if (cached_has_bits & 0x00000008u) {
      GOOGLE_DCHECK(_impl_.a_ != nullptr);
      _impl_.a_->Clear();
    }
    if (cached_has_bits & 0x00000010u) {
      GOOGLE_DCHECK(_impl_.kappa_ != nullptr);
      _impl_.kappa_->Clear();
    }
    if (cached_has_bits & 0x00000020u) {
      GOOGLE_DCHECK(_impl_.dkappa_ != nullptr);
      _impl_.dkappa_->Clear();
    }
  }
  _impl_._has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* PlanningStats::_InternalParse(const char* ptr, ::_pbi::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  while (!ctx->Done(&ptr)) {
    uint32_t tag;
    ptr = ::_pbi::ReadTag(ptr, &tag);
    switch (tag >> 3) {
      // optional .apollo.planning.StatsGroup total_path_length = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 10)) {
          ptr = ctx->ParseMessage(_internal_mutable_total_path_length(), ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.planning.StatsGroup total_path_time = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 18)) {
          ptr = ctx->ParseMessage(_internal_mutable_total_path_time(), ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.planning.StatsGroup v = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 26)) {
          ptr = ctx->ParseMessage(_internal_mutable_v(), ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.planning.StatsGroup a = 4;
      case 4:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 34)) {
          ptr = ctx->ParseMessage(_internal_mutable_a(), ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.planning.StatsGroup kappa = 5;
      case 5:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 42)) {
          ptr = ctx->ParseMessage(_internal_mutable_kappa(), ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.planning.StatsGroup dkappa = 6;
      case 6:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 50)) {
          ptr = ctx->ParseMessage(_internal_mutable_dkappa(), ptr);
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

uint8_t* PlanningStats::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:apollo.planning.PlanningStats)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  // optional .apollo.planning.StatsGroup total_path_length = 1;
  if (cached_has_bits & 0x00000001u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(1, _Internal::total_path_length(this),
        _Internal::total_path_length(this).GetCachedSize(), target, stream);
  }

  // optional .apollo.planning.StatsGroup total_path_time = 2;
  if (cached_has_bits & 0x00000002u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(2, _Internal::total_path_time(this),
        _Internal::total_path_time(this).GetCachedSize(), target, stream);
  }

  // optional .apollo.planning.StatsGroup v = 3;
  if (cached_has_bits & 0x00000004u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(3, _Internal::v(this),
        _Internal::v(this).GetCachedSize(), target, stream);
  }

  // optional .apollo.planning.StatsGroup a = 4;
  if (cached_has_bits & 0x00000008u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(4, _Internal::a(this),
        _Internal::a(this).GetCachedSize(), target, stream);
  }

  // optional .apollo.planning.StatsGroup kappa = 5;
  if (cached_has_bits & 0x00000010u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(5, _Internal::kappa(this),
        _Internal::kappa(this).GetCachedSize(), target, stream);
  }

  // optional .apollo.planning.StatsGroup dkappa = 6;
  if (cached_has_bits & 0x00000020u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(6, _Internal::dkappa(this),
        _Internal::dkappa(this).GetCachedSize(), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:apollo.planning.PlanningStats)
  return target;
}

size_t PlanningStats::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:apollo.planning.PlanningStats)
  size_t total_size = 0;

  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x0000003fu) {
    // optional .apollo.planning.StatsGroup total_path_length = 1;
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.total_path_length_);
    }

    // optional .apollo.planning.StatsGroup total_path_time = 2;
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.total_path_time_);
    }

    // optional .apollo.planning.StatsGroup v = 3;
    if (cached_has_bits & 0x00000004u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.v_);
    }

    // optional .apollo.planning.StatsGroup a = 4;
    if (cached_has_bits & 0x00000008u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.a_);
    }

    // optional .apollo.planning.StatsGroup kappa = 5;
    if (cached_has_bits & 0x00000010u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.kappa_);
    }

    // optional .apollo.planning.StatsGroup dkappa = 6;
    if (cached_has_bits & 0x00000020u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.dkappa_);
    }

  }
  return MaybeComputeUnknownFieldsSize(total_size, &_impl_._cached_size_);
}

const ::PROTOBUF_NAMESPACE_ID::Message::ClassData PlanningStats::_class_data_ = {
    ::PROTOBUF_NAMESPACE_ID::Message::CopyWithSizeCheck,
    PlanningStats::MergeImpl
};
const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*PlanningStats::GetClassData() const { return &_class_data_; }

void PlanningStats::MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to,
                      const ::PROTOBUF_NAMESPACE_ID::Message& from) {
  static_cast<PlanningStats *>(to)->MergeFrom(
      static_cast<const PlanningStats &>(from));
}


void PlanningStats::MergeFrom(const PlanningStats& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:apollo.planning.PlanningStats)
  GOOGLE_DCHECK_NE(&from, this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._impl_._has_bits_[0];
  if (cached_has_bits & 0x0000003fu) {
    if (cached_has_bits & 0x00000001u) {
      _internal_mutable_total_path_length()->::apollo::planning::StatsGroup::MergeFrom(from._internal_total_path_length());
    }
    if (cached_has_bits & 0x00000002u) {
      _internal_mutable_total_path_time()->::apollo::planning::StatsGroup::MergeFrom(from._internal_total_path_time());
    }
    if (cached_has_bits & 0x00000004u) {
      _internal_mutable_v()->::apollo::planning::StatsGroup::MergeFrom(from._internal_v());
    }
    if (cached_has_bits & 0x00000008u) {
      _internal_mutable_a()->::apollo::planning::StatsGroup::MergeFrom(from._internal_a());
    }
    if (cached_has_bits & 0x00000010u) {
      _internal_mutable_kappa()->::apollo::planning::StatsGroup::MergeFrom(from._internal_kappa());
    }
    if (cached_has_bits & 0x00000020u) {
      _internal_mutable_dkappa()->::apollo::planning::StatsGroup::MergeFrom(from._internal_dkappa());
    }
  }
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void PlanningStats::CopyFrom(const PlanningStats& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:apollo.planning.PlanningStats)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool PlanningStats::IsInitialized() const {
  return true;
}

void PlanningStats::InternalSwap(PlanningStats* other) {
  using std::swap;
  _internal_metadata_.InternalSwap(&other->_internal_metadata_);
  swap(_impl_._has_bits_[0], other->_impl_._has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(PlanningStats, _impl_.dkappa_)
      + sizeof(PlanningStats::_impl_.dkappa_)
      - PROTOBUF_FIELD_OFFSET(PlanningStats, _impl_.total_path_length_)>(
          reinterpret_cast<char*>(&_impl_.total_path_length_),
          reinterpret_cast<char*>(&other->_impl_.total_path_length_));
}

::PROTOBUF_NAMESPACE_ID::Metadata PlanningStats::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto_getter, &descriptor_table_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto_once,
      file_level_metadata_modules_2fplanning_2fproto_2fplanning_5fstats_2eproto[1]);
}

// @@protoc_insertion_point(namespace_scope)
}  // namespace planning
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::apollo::planning::StatsGroup*
Arena::CreateMaybeMessage< ::apollo::planning::StatsGroup >(Arena* arena) {
  return Arena::CreateMessageInternal< ::apollo::planning::StatsGroup >(arena);
}
template<> PROTOBUF_NOINLINE ::apollo::planning::PlanningStats*
Arena::CreateMaybeMessage< ::apollo::planning::PlanningStats >(Arena* arena) {
  return Arena::CreateMessageInternal< ::apollo::planning::PlanningStats >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
