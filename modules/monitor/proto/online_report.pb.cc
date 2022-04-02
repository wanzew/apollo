// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/monitor/proto/online_report.proto

#include "modules/monitor/proto/online_report.pb.h"

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
namespace monitor {
PROTOBUF_CONSTEXPR VehicleStateReport::VehicleStateReport(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_._has_bits_)*/{}
  , /*decltype(_impl_._cached_size_)*/{}
  , /*decltype(_impl_.vehicle_info_)*/nullptr
  , /*decltype(_impl_.vehicle_state_)*/nullptr} {}
struct VehicleStateReportDefaultTypeInternal {
  PROTOBUF_CONSTEXPR VehicleStateReportDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~VehicleStateReportDefaultTypeInternal() {}
  union {
    VehicleStateReport _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 VehicleStateReportDefaultTypeInternal _VehicleStateReport_default_instance_;
}  // namespace monitor
}  // namespace apollo
static ::_pb::Metadata file_level_metadata_modules_2fmonitor_2fproto_2fonline_5freport_2eproto[1];
static constexpr ::_pb::EnumDescriptor const** file_level_enum_descriptors_modules_2fmonitor_2fproto_2fonline_5freport_2eproto = nullptr;
static constexpr ::_pb::ServiceDescriptor const** file_level_service_descriptors_modules_2fmonitor_2fproto_2fonline_5freport_2eproto = nullptr;

const uint32_t TableStruct_modules_2fmonitor_2fproto_2fonline_5freport_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::apollo::monitor::VehicleStateReport, _impl_._has_bits_),
  PROTOBUF_FIELD_OFFSET(::apollo::monitor::VehicleStateReport, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::apollo::monitor::VehicleStateReport, _impl_.vehicle_info_),
  PROTOBUF_FIELD_OFFSET(::apollo::monitor::VehicleStateReport, _impl_.vehicle_state_),
  0,
  1,
};
static const ::_pbi::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 8, -1, sizeof(::apollo::monitor::VehicleStateReport)},
};

static const ::_pb::Message* const file_default_instances[] = {
  &::apollo::monitor::_VehicleStateReport_default_instance_._instance,
};

const char descriptor_table_protodef_modules_2fmonitor_2fproto_2fonline_5freport_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n)modules/monitor/proto/online_report.pr"
  "oto\022\016apollo.monitor\0326modules/common/vehi"
  "cle_state/proto/vehicle_state.proto\032$mod"
  "ules/data/proto/static_info.proto\"x\n\022Veh"
  "icleStateReport\022.\n\014vehicle_info\030\001 \001(\0132\030."
  "apollo.data.VehicleInfo\0222\n\rvehicle_state"
  "\030\002 \001(\0132\033.apollo.common.VehicleState"
  ;
static const ::_pbi::DescriptorTable* const descriptor_table_modules_2fmonitor_2fproto_2fonline_5freport_2eproto_deps[2] = {
  &::descriptor_table_modules_2fcommon_2fvehicle_5fstate_2fproto_2fvehicle_5fstate_2eproto,
  &::descriptor_table_modules_2fdata_2fproto_2fstatic_5finfo_2eproto,
};
static ::_pbi::once_flag descriptor_table_modules_2fmonitor_2fproto_2fonline_5freport_2eproto_once;
const ::_pbi::DescriptorTable descriptor_table_modules_2fmonitor_2fproto_2fonline_5freport_2eproto = {
    false, false, 275, descriptor_table_protodef_modules_2fmonitor_2fproto_2fonline_5freport_2eproto,
    "modules/monitor/proto/online_report.proto",
    &descriptor_table_modules_2fmonitor_2fproto_2fonline_5freport_2eproto_once, descriptor_table_modules_2fmonitor_2fproto_2fonline_5freport_2eproto_deps, 2, 1,
    schemas, file_default_instances, TableStruct_modules_2fmonitor_2fproto_2fonline_5freport_2eproto::offsets,
    file_level_metadata_modules_2fmonitor_2fproto_2fonline_5freport_2eproto, file_level_enum_descriptors_modules_2fmonitor_2fproto_2fonline_5freport_2eproto,
    file_level_service_descriptors_modules_2fmonitor_2fproto_2fonline_5freport_2eproto,
};
PROTOBUF_ATTRIBUTE_WEAK const ::_pbi::DescriptorTable* descriptor_table_modules_2fmonitor_2fproto_2fonline_5freport_2eproto_getter() {
  return &descriptor_table_modules_2fmonitor_2fproto_2fonline_5freport_2eproto;
}

// Force running AddDescriptors() at dynamic initialization time.
PROTOBUF_ATTRIBUTE_INIT_PRIORITY2 static ::_pbi::AddDescriptorsRunner dynamic_init_dummy_modules_2fmonitor_2fproto_2fonline_5freport_2eproto(&descriptor_table_modules_2fmonitor_2fproto_2fonline_5freport_2eproto);
namespace apollo {
namespace monitor {

// ===================================================================

class VehicleStateReport::_Internal {
 public:
  using HasBits = decltype(std::declval<VehicleStateReport>()._impl_._has_bits_);
  static const ::apollo::data::VehicleInfo& vehicle_info(const VehicleStateReport* msg);
  static void set_has_vehicle_info(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static const ::apollo::common::VehicleState& vehicle_state(const VehicleStateReport* msg);
  static void set_has_vehicle_state(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
};

const ::apollo::data::VehicleInfo&
VehicleStateReport::_Internal::vehicle_info(const VehicleStateReport* msg) {
  return *msg->_impl_.vehicle_info_;
}
const ::apollo::common::VehicleState&
VehicleStateReport::_Internal::vehicle_state(const VehicleStateReport* msg) {
  return *msg->_impl_.vehicle_state_;
}
void VehicleStateReport::clear_vehicle_info() {
  if (_impl_.vehicle_info_ != nullptr) _impl_.vehicle_info_->Clear();
  _impl_._has_bits_[0] &= ~0x00000001u;
}
void VehicleStateReport::clear_vehicle_state() {
  if (_impl_.vehicle_state_ != nullptr) _impl_.vehicle_state_->Clear();
  _impl_._has_bits_[0] &= ~0x00000002u;
}
VehicleStateReport::VehicleStateReport(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:apollo.monitor.VehicleStateReport)
}
VehicleStateReport::VehicleStateReport(const VehicleStateReport& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){from._impl_._has_bits_}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.vehicle_info_){nullptr}
    , decltype(_impl_.vehicle_state_){nullptr}};

  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  if (from._internal_has_vehicle_info()) {
    _impl_.vehicle_info_ = new ::apollo::data::VehicleInfo(*from._impl_.vehicle_info_);
  }
  if (from._internal_has_vehicle_state()) {
    _impl_.vehicle_state_ = new ::apollo::common::VehicleState(*from._impl_.vehicle_state_);
  }
  // @@protoc_insertion_point(copy_constructor:apollo.monitor.VehicleStateReport)
}

inline void VehicleStateReport::SharedCtor(
    ::_pb::Arena* arena, bool is_message_owned) {
  (void)arena;
  (void)is_message_owned;
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.vehicle_info_){nullptr}
    , decltype(_impl_.vehicle_state_){nullptr}
  };
}

VehicleStateReport::~VehicleStateReport() {
  // @@protoc_insertion_point(destructor:apollo.monitor.VehicleStateReport)
  if (auto *arena = _internal_metadata_.DeleteReturnArena<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>()) {
  (void)arena;
    return;
  }
  SharedDtor();
}

inline void VehicleStateReport::SharedDtor() {
  GOOGLE_DCHECK(GetArenaForAllocation() == nullptr);
  if (this != internal_default_instance()) delete _impl_.vehicle_info_;
  if (this != internal_default_instance()) delete _impl_.vehicle_state_;
}

void VehicleStateReport::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void VehicleStateReport::Clear() {
// @@protoc_insertion_point(message_clear_start:apollo.monitor.VehicleStateReport)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    if (cached_has_bits & 0x00000001u) {
      GOOGLE_DCHECK(_impl_.vehicle_info_ != nullptr);
      _impl_.vehicle_info_->Clear();
    }
    if (cached_has_bits & 0x00000002u) {
      GOOGLE_DCHECK(_impl_.vehicle_state_ != nullptr);
      _impl_.vehicle_state_->Clear();
    }
  }
  _impl_._has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* VehicleStateReport::_InternalParse(const char* ptr, ::_pbi::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  while (!ctx->Done(&ptr)) {
    uint32_t tag;
    ptr = ::_pbi::ReadTag(ptr, &tag);
    switch (tag >> 3) {
      // optional .apollo.data.VehicleInfo vehicle_info = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 10)) {
          ptr = ctx->ParseMessage(_internal_mutable_vehicle_info(), ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.common.VehicleState vehicle_state = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 18)) {
          ptr = ctx->ParseMessage(_internal_mutable_vehicle_state(), ptr);
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

uint8_t* VehicleStateReport::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:apollo.monitor.VehicleStateReport)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  // optional .apollo.data.VehicleInfo vehicle_info = 1;
  if (cached_has_bits & 0x00000001u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(1, _Internal::vehicle_info(this),
        _Internal::vehicle_info(this).GetCachedSize(), target, stream);
  }

  // optional .apollo.common.VehicleState vehicle_state = 2;
  if (cached_has_bits & 0x00000002u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(2, _Internal::vehicle_state(this),
        _Internal::vehicle_state(this).GetCachedSize(), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:apollo.monitor.VehicleStateReport)
  return target;
}

size_t VehicleStateReport::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:apollo.monitor.VehicleStateReport)
  size_t total_size = 0;

  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    // optional .apollo.data.VehicleInfo vehicle_info = 1;
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.vehicle_info_);
    }

    // optional .apollo.common.VehicleState vehicle_state = 2;
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.vehicle_state_);
    }

  }
  return MaybeComputeUnknownFieldsSize(total_size, &_impl_._cached_size_);
}

const ::PROTOBUF_NAMESPACE_ID::Message::ClassData VehicleStateReport::_class_data_ = {
    ::PROTOBUF_NAMESPACE_ID::Message::CopyWithSizeCheck,
    VehicleStateReport::MergeImpl
};
const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*VehicleStateReport::GetClassData() const { return &_class_data_; }

void VehicleStateReport::MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to,
                      const ::PROTOBUF_NAMESPACE_ID::Message& from) {
  static_cast<VehicleStateReport *>(to)->MergeFrom(
      static_cast<const VehicleStateReport &>(from));
}


void VehicleStateReport::MergeFrom(const VehicleStateReport& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:apollo.monitor.VehicleStateReport)
  GOOGLE_DCHECK_NE(&from, this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._impl_._has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    if (cached_has_bits & 0x00000001u) {
      _internal_mutable_vehicle_info()->::apollo::data::VehicleInfo::MergeFrom(from._internal_vehicle_info());
    }
    if (cached_has_bits & 0x00000002u) {
      _internal_mutable_vehicle_state()->::apollo::common::VehicleState::MergeFrom(from._internal_vehicle_state());
    }
  }
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void VehicleStateReport::CopyFrom(const VehicleStateReport& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:apollo.monitor.VehicleStateReport)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool VehicleStateReport::IsInitialized() const {
  return true;
}

void VehicleStateReport::InternalSwap(VehicleStateReport* other) {
  using std::swap;
  _internal_metadata_.InternalSwap(&other->_internal_metadata_);
  swap(_impl_._has_bits_[0], other->_impl_._has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(VehicleStateReport, _impl_.vehicle_state_)
      + sizeof(VehicleStateReport::_impl_.vehicle_state_)
      - PROTOBUF_FIELD_OFFSET(VehicleStateReport, _impl_.vehicle_info_)>(
          reinterpret_cast<char*>(&_impl_.vehicle_info_),
          reinterpret_cast<char*>(&other->_impl_.vehicle_info_));
}

::PROTOBUF_NAMESPACE_ID::Metadata VehicleStateReport::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_modules_2fmonitor_2fproto_2fonline_5freport_2eproto_getter, &descriptor_table_modules_2fmonitor_2fproto_2fonline_5freport_2eproto_once,
      file_level_metadata_modules_2fmonitor_2fproto_2fonline_5freport_2eproto[0]);
}

// @@protoc_insertion_point(namespace_scope)
}  // namespace monitor
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::apollo::monitor::VehicleStateReport*
Arena::CreateMaybeMessage< ::apollo::monitor::VehicleStateReport >(Arena* arena) {
  return Arena::CreateMessageInternal< ::apollo::monitor::VehicleStateReport >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
