// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/control/proto/input_debug.proto

#include "modules/control/proto/input_debug.pb.h"

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
namespace control {
PROTOBUF_CONSTEXPR InputDebug::InputDebug(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_._has_bits_)*/{}
  , /*decltype(_impl_._cached_size_)*/{}
  , /*decltype(_impl_.localization_header_)*/nullptr
  , /*decltype(_impl_.canbus_header_)*/nullptr
  , /*decltype(_impl_.trajectory_header_)*/nullptr
  , /*decltype(_impl_.latest_replan_trajectory_header_)*/nullptr} {}
struct InputDebugDefaultTypeInternal {
  PROTOBUF_CONSTEXPR InputDebugDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~InputDebugDefaultTypeInternal() {}
  union {
    InputDebug _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 InputDebugDefaultTypeInternal _InputDebug_default_instance_;
}  // namespace control
}  // namespace apollo
static ::_pb::Metadata file_level_metadata_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto[1];
static constexpr ::_pb::EnumDescriptor const** file_level_enum_descriptors_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto = nullptr;
static constexpr ::_pb::ServiceDescriptor const** file_level_service_descriptors_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto = nullptr;

const uint32_t TableStruct_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::apollo::control::InputDebug, _impl_._has_bits_),
  PROTOBUF_FIELD_OFFSET(::apollo::control::InputDebug, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::apollo::control::InputDebug, _impl_.localization_header_),
  PROTOBUF_FIELD_OFFSET(::apollo::control::InputDebug, _impl_.canbus_header_),
  PROTOBUF_FIELD_OFFSET(::apollo::control::InputDebug, _impl_.trajectory_header_),
  PROTOBUF_FIELD_OFFSET(::apollo::control::InputDebug, _impl_.latest_replan_trajectory_header_),
  0,
  1,
  2,
  3,
};
static const ::_pbi::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 10, -1, sizeof(::apollo::control::InputDebug)},
};

static const ::_pb::Message* const file_default_instances[] = {
  &::apollo::control::_InputDebug_default_instance_._instance,
};

const char descriptor_table_protodef_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n\'modules/control/proto/input_debug.prot"
  "o\022\016apollo.control\032!modules/common/proto/"
  "header.proto\"\340\001\n\nInputDebug\0222\n\023localizat"
  "ion_header\030\001 \001(\0132\025.apollo.common.Header\022"
  ",\n\rcanbus_header\030\002 \001(\0132\025.apollo.common.H"
  "eader\0220\n\021trajectory_header\030\003 \001(\0132\025.apoll"
  "o.common.Header\022>\n\037latest_replan_traject"
  "ory_header\030\004 \001(\0132\025.apollo.common.Header"
  ;
static const ::_pbi::DescriptorTable* const descriptor_table_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto_deps[1] = {
  &::descriptor_table_modules_2fcommon_2fproto_2fheader_2eproto,
};
static ::_pbi::once_flag descriptor_table_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto_once;
const ::_pbi::DescriptorTable descriptor_table_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto = {
    false, false, 319, descriptor_table_protodef_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto,
    "modules/control/proto/input_debug.proto",
    &descriptor_table_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto_once, descriptor_table_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto_deps, 1, 1,
    schemas, file_default_instances, TableStruct_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto::offsets,
    file_level_metadata_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto, file_level_enum_descriptors_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto,
    file_level_service_descriptors_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto,
};
PROTOBUF_ATTRIBUTE_WEAK const ::_pbi::DescriptorTable* descriptor_table_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto_getter() {
  return &descriptor_table_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto;
}

// Force running AddDescriptors() at dynamic initialization time.
PROTOBUF_ATTRIBUTE_INIT_PRIORITY2 static ::_pbi::AddDescriptorsRunner dynamic_init_dummy_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto(&descriptor_table_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto);
namespace apollo {
namespace control {

// ===================================================================

class InputDebug::_Internal {
 public:
  using HasBits = decltype(std::declval<InputDebug>()._impl_._has_bits_);
  static const ::apollo::common::Header& localization_header(const InputDebug* msg);
  static void set_has_localization_header(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static const ::apollo::common::Header& canbus_header(const InputDebug* msg);
  static void set_has_canbus_header(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static const ::apollo::common::Header& trajectory_header(const InputDebug* msg);
  static void set_has_trajectory_header(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static const ::apollo::common::Header& latest_replan_trajectory_header(const InputDebug* msg);
  static void set_has_latest_replan_trajectory_header(HasBits* has_bits) {
    (*has_bits)[0] |= 8u;
  }
};

const ::apollo::common::Header&
InputDebug::_Internal::localization_header(const InputDebug* msg) {
  return *msg->_impl_.localization_header_;
}
const ::apollo::common::Header&
InputDebug::_Internal::canbus_header(const InputDebug* msg) {
  return *msg->_impl_.canbus_header_;
}
const ::apollo::common::Header&
InputDebug::_Internal::trajectory_header(const InputDebug* msg) {
  return *msg->_impl_.trajectory_header_;
}
const ::apollo::common::Header&
InputDebug::_Internal::latest_replan_trajectory_header(const InputDebug* msg) {
  return *msg->_impl_.latest_replan_trajectory_header_;
}
void InputDebug::clear_localization_header() {
  if (_impl_.localization_header_ != nullptr) _impl_.localization_header_->Clear();
  _impl_._has_bits_[0] &= ~0x00000001u;
}
void InputDebug::clear_canbus_header() {
  if (_impl_.canbus_header_ != nullptr) _impl_.canbus_header_->Clear();
  _impl_._has_bits_[0] &= ~0x00000002u;
}
void InputDebug::clear_trajectory_header() {
  if (_impl_.trajectory_header_ != nullptr) _impl_.trajectory_header_->Clear();
  _impl_._has_bits_[0] &= ~0x00000004u;
}
void InputDebug::clear_latest_replan_trajectory_header() {
  if (_impl_.latest_replan_trajectory_header_ != nullptr) _impl_.latest_replan_trajectory_header_->Clear();
  _impl_._has_bits_[0] &= ~0x00000008u;
}
InputDebug::InputDebug(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:apollo.control.InputDebug)
}
InputDebug::InputDebug(const InputDebug& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){from._impl_._has_bits_}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.localization_header_){nullptr}
    , decltype(_impl_.canbus_header_){nullptr}
    , decltype(_impl_.trajectory_header_){nullptr}
    , decltype(_impl_.latest_replan_trajectory_header_){nullptr}};

  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  if (from._internal_has_localization_header()) {
    _impl_.localization_header_ = new ::apollo::common::Header(*from._impl_.localization_header_);
  }
  if (from._internal_has_canbus_header()) {
    _impl_.canbus_header_ = new ::apollo::common::Header(*from._impl_.canbus_header_);
  }
  if (from._internal_has_trajectory_header()) {
    _impl_.trajectory_header_ = new ::apollo::common::Header(*from._impl_.trajectory_header_);
  }
  if (from._internal_has_latest_replan_trajectory_header()) {
    _impl_.latest_replan_trajectory_header_ = new ::apollo::common::Header(*from._impl_.latest_replan_trajectory_header_);
  }
  // @@protoc_insertion_point(copy_constructor:apollo.control.InputDebug)
}

inline void InputDebug::SharedCtor(
    ::_pb::Arena* arena, bool is_message_owned) {
  (void)arena;
  (void)is_message_owned;
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.localization_header_){nullptr}
    , decltype(_impl_.canbus_header_){nullptr}
    , decltype(_impl_.trajectory_header_){nullptr}
    , decltype(_impl_.latest_replan_trajectory_header_){nullptr}
  };
}

InputDebug::~InputDebug() {
  // @@protoc_insertion_point(destructor:apollo.control.InputDebug)
  if (auto *arena = _internal_metadata_.DeleteReturnArena<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>()) {
  (void)arena;
    return;
  }
  SharedDtor();
}

inline void InputDebug::SharedDtor() {
  GOOGLE_DCHECK(GetArenaForAllocation() == nullptr);
  if (this != internal_default_instance()) delete _impl_.localization_header_;
  if (this != internal_default_instance()) delete _impl_.canbus_header_;
  if (this != internal_default_instance()) delete _impl_.trajectory_header_;
  if (this != internal_default_instance()) delete _impl_.latest_replan_trajectory_header_;
}

void InputDebug::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void InputDebug::Clear() {
// @@protoc_insertion_point(message_clear_start:apollo.control.InputDebug)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x0000000fu) {
    if (cached_has_bits & 0x00000001u) {
      GOOGLE_DCHECK(_impl_.localization_header_ != nullptr);
      _impl_.localization_header_->Clear();
    }
    if (cached_has_bits & 0x00000002u) {
      GOOGLE_DCHECK(_impl_.canbus_header_ != nullptr);
      _impl_.canbus_header_->Clear();
    }
    if (cached_has_bits & 0x00000004u) {
      GOOGLE_DCHECK(_impl_.trajectory_header_ != nullptr);
      _impl_.trajectory_header_->Clear();
    }
    if (cached_has_bits & 0x00000008u) {
      GOOGLE_DCHECK(_impl_.latest_replan_trajectory_header_ != nullptr);
      _impl_.latest_replan_trajectory_header_->Clear();
    }
  }
  _impl_._has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* InputDebug::_InternalParse(const char* ptr, ::_pbi::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  while (!ctx->Done(&ptr)) {
    uint32_t tag;
    ptr = ::_pbi::ReadTag(ptr, &tag);
    switch (tag >> 3) {
      // optional .apollo.common.Header localization_header = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 10)) {
          ptr = ctx->ParseMessage(_internal_mutable_localization_header(), ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.common.Header canbus_header = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 18)) {
          ptr = ctx->ParseMessage(_internal_mutable_canbus_header(), ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.common.Header trajectory_header = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 26)) {
          ptr = ctx->ParseMessage(_internal_mutable_trajectory_header(), ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.common.Header latest_replan_trajectory_header = 4;
      case 4:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 34)) {
          ptr = ctx->ParseMessage(_internal_mutable_latest_replan_trajectory_header(), ptr);
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

uint8_t* InputDebug::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:apollo.control.InputDebug)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  // optional .apollo.common.Header localization_header = 1;
  if (cached_has_bits & 0x00000001u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(1, _Internal::localization_header(this),
        _Internal::localization_header(this).GetCachedSize(), target, stream);
  }

  // optional .apollo.common.Header canbus_header = 2;
  if (cached_has_bits & 0x00000002u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(2, _Internal::canbus_header(this),
        _Internal::canbus_header(this).GetCachedSize(), target, stream);
  }

  // optional .apollo.common.Header trajectory_header = 3;
  if (cached_has_bits & 0x00000004u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(3, _Internal::trajectory_header(this),
        _Internal::trajectory_header(this).GetCachedSize(), target, stream);
  }

  // optional .apollo.common.Header latest_replan_trajectory_header = 4;
  if (cached_has_bits & 0x00000008u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(4, _Internal::latest_replan_trajectory_header(this),
        _Internal::latest_replan_trajectory_header(this).GetCachedSize(), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:apollo.control.InputDebug)
  return target;
}

size_t InputDebug::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:apollo.control.InputDebug)
  size_t total_size = 0;

  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x0000000fu) {
    // optional .apollo.common.Header localization_header = 1;
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.localization_header_);
    }

    // optional .apollo.common.Header canbus_header = 2;
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.canbus_header_);
    }

    // optional .apollo.common.Header trajectory_header = 3;
    if (cached_has_bits & 0x00000004u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.trajectory_header_);
    }

    // optional .apollo.common.Header latest_replan_trajectory_header = 4;
    if (cached_has_bits & 0x00000008u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.latest_replan_trajectory_header_);
    }

  }
  return MaybeComputeUnknownFieldsSize(total_size, &_impl_._cached_size_);
}

const ::PROTOBUF_NAMESPACE_ID::Message::ClassData InputDebug::_class_data_ = {
    ::PROTOBUF_NAMESPACE_ID::Message::CopyWithSizeCheck,
    InputDebug::MergeImpl
};
const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*InputDebug::GetClassData() const { return &_class_data_; }

void InputDebug::MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to,
                      const ::PROTOBUF_NAMESPACE_ID::Message& from) {
  static_cast<InputDebug *>(to)->MergeFrom(
      static_cast<const InputDebug &>(from));
}


void InputDebug::MergeFrom(const InputDebug& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:apollo.control.InputDebug)
  GOOGLE_DCHECK_NE(&from, this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._impl_._has_bits_[0];
  if (cached_has_bits & 0x0000000fu) {
    if (cached_has_bits & 0x00000001u) {
      _internal_mutable_localization_header()->::apollo::common::Header::MergeFrom(from._internal_localization_header());
    }
    if (cached_has_bits & 0x00000002u) {
      _internal_mutable_canbus_header()->::apollo::common::Header::MergeFrom(from._internal_canbus_header());
    }
    if (cached_has_bits & 0x00000004u) {
      _internal_mutable_trajectory_header()->::apollo::common::Header::MergeFrom(from._internal_trajectory_header());
    }
    if (cached_has_bits & 0x00000008u) {
      _internal_mutable_latest_replan_trajectory_header()->::apollo::common::Header::MergeFrom(from._internal_latest_replan_trajectory_header());
    }
  }
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void InputDebug::CopyFrom(const InputDebug& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:apollo.control.InputDebug)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool InputDebug::IsInitialized() const {
  return true;
}

void InputDebug::InternalSwap(InputDebug* other) {
  using std::swap;
  _internal_metadata_.InternalSwap(&other->_internal_metadata_);
  swap(_impl_._has_bits_[0], other->_impl_._has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(InputDebug, _impl_.latest_replan_trajectory_header_)
      + sizeof(InputDebug::_impl_.latest_replan_trajectory_header_)
      - PROTOBUF_FIELD_OFFSET(InputDebug, _impl_.localization_header_)>(
          reinterpret_cast<char*>(&_impl_.localization_header_),
          reinterpret_cast<char*>(&other->_impl_.localization_header_));
}

::PROTOBUF_NAMESPACE_ID::Metadata InputDebug::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto_getter, &descriptor_table_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto_once,
      file_level_metadata_modules_2fcontrol_2fproto_2finput_5fdebug_2eproto[0]);
}

// @@protoc_insertion_point(namespace_scope)
}  // namespace control
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::apollo::control::InputDebug*
Arena::CreateMaybeMessage< ::apollo::control::InputDebug >(Arena* arena) {
  return Arena::CreateMessageInternal< ::apollo::control::InputDebug >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
