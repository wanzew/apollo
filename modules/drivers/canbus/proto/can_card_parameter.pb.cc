// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/drivers/canbus/proto/can_card_parameter.proto

#include "modules/drivers/canbus/proto/can_card_parameter.pb.h"

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
namespace drivers {
namespace canbus {
PROTOBUF_CONSTEXPR CANCardParameter::CANCardParameter(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_._has_bits_)*/{}
  , /*decltype(_impl_._cached_size_)*/{}
  , /*decltype(_impl_.brand_)*/0
  , /*decltype(_impl_.type_)*/0
  , /*decltype(_impl_.channel_id_)*/0
  , /*decltype(_impl_.interface_)*/0
  , /*decltype(_impl_.num_ports_)*/4u} {}
struct CANCardParameterDefaultTypeInternal {
  PROTOBUF_CONSTEXPR CANCardParameterDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~CANCardParameterDefaultTypeInternal() {}
  union {
    CANCardParameter _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 CANCardParameterDefaultTypeInternal _CANCardParameter_default_instance_;
}  // namespace canbus
}  // namespace drivers
}  // namespace apollo
static ::_pb::Metadata file_level_metadata_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto[1];
static const ::_pb::EnumDescriptor* file_level_enum_descriptors_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto[4];
static constexpr ::_pb::ServiceDescriptor const** file_level_service_descriptors_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto = nullptr;

const uint32_t TableStruct_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::apollo::drivers::canbus::CANCardParameter, _impl_._has_bits_),
  PROTOBUF_FIELD_OFFSET(::apollo::drivers::canbus::CANCardParameter, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::apollo::drivers::canbus::CANCardParameter, _impl_.brand_),
  PROTOBUF_FIELD_OFFSET(::apollo::drivers::canbus::CANCardParameter, _impl_.type_),
  PROTOBUF_FIELD_OFFSET(::apollo::drivers::canbus::CANCardParameter, _impl_.channel_id_),
  PROTOBUF_FIELD_OFFSET(::apollo::drivers::canbus::CANCardParameter, _impl_.interface_),
  PROTOBUF_FIELD_OFFSET(::apollo::drivers::canbus::CANCardParameter, _impl_.num_ports_),
  0,
  1,
  2,
  3,
  4,
};
static const ::_pbi::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 11, -1, sizeof(::apollo::drivers::canbus::CANCardParameter)},
};

static const ::_pb::Message* const file_default_instances[] = {
  &::apollo::drivers::canbus::_CANCardParameter_default_instance_._instance,
};

const char descriptor_table_protodef_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n5modules/drivers/canbus/proto/can_card_"
  "parameter.proto\022\025apollo.drivers.canbus\"\251"
  "\005\n\020CANCardParameter\022C\n\005brand\030\001 \001(\01624.apo"
  "llo.drivers.canbus.CANCardParameter.CANC"
  "ardBrand\022A\n\004type\030\002 \001(\01623.apollo.drivers."
  "canbus.CANCardParameter.CANCardType\022H\n\nc"
  "hannel_id\030\003 \001(\01624.apollo.drivers.canbus."
  "CANCardParameter.CANChannelId\022G\n\tinterfa"
  "ce\030\004 \001(\01624.apollo.drivers.canbus.CANCard"
  "Parameter.CANInterface\022\024\n\tnum_ports\030\005 \001("
  "\r:\0014\"M\n\014CANCardBrand\022\014\n\010FAKE_CAN\020\000\022\013\n\007ES"
  "D_CAN\020\001\022\022\n\016SOCKET_CAN_RAW\020\002\022\016\n\nHERMES_CA"
  "N\020\003\")\n\013CANCardType\022\014\n\010PCI_CARD\020\000\022\014\n\010USB_"
  "CARD\020\001\"\265\001\n\014CANChannelId\022\023\n\017CHANNEL_ID_ZE"
  "RO\020\000\022\022\n\016CHANNEL_ID_ONE\020\001\022\022\n\016CHANNEL_ID_T"
  "WO\020\002\022\024\n\020CHANNEL_ID_THREE\020\003\022\023\n\017CHANNEL_ID"
  "_FOUR\020\004\022\023\n\017CHANNEL_ID_FIVE\020\005\022\022\n\016CHANNEL_"
  "ID_SIX\020\006\022\024\n\020CHANNEL_ID_SEVEN\020\007\"2\n\014CANInt"
  "erface\022\n\n\006NATIVE\020\000\022\013\n\007VIRTUAL\020\001\022\t\n\005SLCAN"
  "\020\002"
  ;
static ::_pbi::once_flag descriptor_table_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto_once;
const ::_pbi::DescriptorTable descriptor_table_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto = {
    false, false, 762, descriptor_table_protodef_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto,
    "modules/drivers/canbus/proto/can_card_parameter.proto",
    &descriptor_table_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto_once, nullptr, 0, 1,
    schemas, file_default_instances, TableStruct_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto::offsets,
    file_level_metadata_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto, file_level_enum_descriptors_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto,
    file_level_service_descriptors_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto,
};
PROTOBUF_ATTRIBUTE_WEAK const ::_pbi::DescriptorTable* descriptor_table_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto_getter() {
  return &descriptor_table_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto;
}

// Force running AddDescriptors() at dynamic initialization time.
PROTOBUF_ATTRIBUTE_INIT_PRIORITY2 static ::_pbi::AddDescriptorsRunner dynamic_init_dummy_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto(&descriptor_table_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto);
namespace apollo {
namespace drivers {
namespace canbus {
const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* CANCardParameter_CANCardBrand_descriptor() {
  ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&descriptor_table_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto);
  return file_level_enum_descriptors_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto[0];
}
bool CANCardParameter_CANCardBrand_IsValid(int value) {
  switch (value) {
    case 0:
    case 1:
    case 2:
    case 3:
      return true;
    default:
      return false;
  }
}

#if (__cplusplus < 201703) && (!defined(_MSC_VER) || (_MSC_VER >= 1900 && _MSC_VER < 1912))
constexpr CANCardParameter_CANCardBrand CANCardParameter::FAKE_CAN;
constexpr CANCardParameter_CANCardBrand CANCardParameter::ESD_CAN;
constexpr CANCardParameter_CANCardBrand CANCardParameter::SOCKET_CAN_RAW;
constexpr CANCardParameter_CANCardBrand CANCardParameter::HERMES_CAN;
constexpr CANCardParameter_CANCardBrand CANCardParameter::CANCardBrand_MIN;
constexpr CANCardParameter_CANCardBrand CANCardParameter::CANCardBrand_MAX;
constexpr int CANCardParameter::CANCardBrand_ARRAYSIZE;
#endif  // (__cplusplus < 201703) && (!defined(_MSC_VER) || (_MSC_VER >= 1900 && _MSC_VER < 1912))
const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* CANCardParameter_CANCardType_descriptor() {
  ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&descriptor_table_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto);
  return file_level_enum_descriptors_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto[1];
}
bool CANCardParameter_CANCardType_IsValid(int value) {
  switch (value) {
    case 0:
    case 1:
      return true;
    default:
      return false;
  }
}

#if (__cplusplus < 201703) && (!defined(_MSC_VER) || (_MSC_VER >= 1900 && _MSC_VER < 1912))
constexpr CANCardParameter_CANCardType CANCardParameter::PCI_CARD;
constexpr CANCardParameter_CANCardType CANCardParameter::USB_CARD;
constexpr CANCardParameter_CANCardType CANCardParameter::CANCardType_MIN;
constexpr CANCardParameter_CANCardType CANCardParameter::CANCardType_MAX;
constexpr int CANCardParameter::CANCardType_ARRAYSIZE;
#endif  // (__cplusplus < 201703) && (!defined(_MSC_VER) || (_MSC_VER >= 1900 && _MSC_VER < 1912))
const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* CANCardParameter_CANChannelId_descriptor() {
  ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&descriptor_table_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto);
  return file_level_enum_descriptors_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto[2];
}
bool CANCardParameter_CANChannelId_IsValid(int value) {
  switch (value) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
      return true;
    default:
      return false;
  }
}

#if (__cplusplus < 201703) && (!defined(_MSC_VER) || (_MSC_VER >= 1900 && _MSC_VER < 1912))
constexpr CANCardParameter_CANChannelId CANCardParameter::CHANNEL_ID_ZERO;
constexpr CANCardParameter_CANChannelId CANCardParameter::CHANNEL_ID_ONE;
constexpr CANCardParameter_CANChannelId CANCardParameter::CHANNEL_ID_TWO;
constexpr CANCardParameter_CANChannelId CANCardParameter::CHANNEL_ID_THREE;
constexpr CANCardParameter_CANChannelId CANCardParameter::CHANNEL_ID_FOUR;
constexpr CANCardParameter_CANChannelId CANCardParameter::CHANNEL_ID_FIVE;
constexpr CANCardParameter_CANChannelId CANCardParameter::CHANNEL_ID_SIX;
constexpr CANCardParameter_CANChannelId CANCardParameter::CHANNEL_ID_SEVEN;
constexpr CANCardParameter_CANChannelId CANCardParameter::CANChannelId_MIN;
constexpr CANCardParameter_CANChannelId CANCardParameter::CANChannelId_MAX;
constexpr int CANCardParameter::CANChannelId_ARRAYSIZE;
#endif  // (__cplusplus < 201703) && (!defined(_MSC_VER) || (_MSC_VER >= 1900 && _MSC_VER < 1912))
const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* CANCardParameter_CANInterface_descriptor() {
  ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&descriptor_table_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto);
  return file_level_enum_descriptors_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto[3];
}
bool CANCardParameter_CANInterface_IsValid(int value) {
  switch (value) {
    case 0:
    case 1:
    case 2:
      return true;
    default:
      return false;
  }
}

#if (__cplusplus < 201703) && (!defined(_MSC_VER) || (_MSC_VER >= 1900 && _MSC_VER < 1912))
constexpr CANCardParameter_CANInterface CANCardParameter::NATIVE;
constexpr CANCardParameter_CANInterface CANCardParameter::VIRTUAL;
constexpr CANCardParameter_CANInterface CANCardParameter::SLCAN;
constexpr CANCardParameter_CANInterface CANCardParameter::CANInterface_MIN;
constexpr CANCardParameter_CANInterface CANCardParameter::CANInterface_MAX;
constexpr int CANCardParameter::CANInterface_ARRAYSIZE;
#endif  // (__cplusplus < 201703) && (!defined(_MSC_VER) || (_MSC_VER >= 1900 && _MSC_VER < 1912))

// ===================================================================

class CANCardParameter::_Internal {
 public:
  using HasBits = decltype(std::declval<CANCardParameter>()._impl_._has_bits_);
  static void set_has_brand(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_type(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static void set_has_channel_id(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static void set_has_interface(HasBits* has_bits) {
    (*has_bits)[0] |= 8u;
  }
  static void set_has_num_ports(HasBits* has_bits) {
    (*has_bits)[0] |= 16u;
  }
};

CANCardParameter::CANCardParameter(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:apollo.drivers.canbus.CANCardParameter)
}
CANCardParameter::CANCardParameter(const CANCardParameter& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){from._impl_._has_bits_}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.brand_){}
    , decltype(_impl_.type_){}
    , decltype(_impl_.channel_id_){}
    , decltype(_impl_.interface_){}
    , decltype(_impl_.num_ports_){}};

  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::memcpy(&_impl_.brand_, &from._impl_.brand_,
    static_cast<size_t>(reinterpret_cast<char*>(&_impl_.num_ports_) -
    reinterpret_cast<char*>(&_impl_.brand_)) + sizeof(_impl_.num_ports_));
  // @@protoc_insertion_point(copy_constructor:apollo.drivers.canbus.CANCardParameter)
}

inline void CANCardParameter::SharedCtor(
    ::_pb::Arena* arena, bool is_message_owned) {
  (void)arena;
  (void)is_message_owned;
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.brand_){0}
    , decltype(_impl_.type_){0}
    , decltype(_impl_.channel_id_){0}
    , decltype(_impl_.interface_){0}
    , decltype(_impl_.num_ports_){4u}
  };
}

CANCardParameter::~CANCardParameter() {
  // @@protoc_insertion_point(destructor:apollo.drivers.canbus.CANCardParameter)
  if (auto *arena = _internal_metadata_.DeleteReturnArena<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>()) {
  (void)arena;
    return;
  }
  SharedDtor();
}

inline void CANCardParameter::SharedDtor() {
  GOOGLE_DCHECK(GetArenaForAllocation() == nullptr);
}

void CANCardParameter::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void CANCardParameter::Clear() {
// @@protoc_insertion_point(message_clear_start:apollo.drivers.canbus.CANCardParameter)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    ::memset(&_impl_.brand_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&_impl_.interface_) -
        reinterpret_cast<char*>(&_impl_.brand_)) + sizeof(_impl_.interface_));
    _impl_.num_ports_ = 4u;
  }
  _impl_._has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* CANCardParameter::_InternalParse(const char* ptr, ::_pbi::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  while (!ctx->Done(&ptr)) {
    uint32_t tag;
    ptr = ::_pbi::ReadTag(ptr, &tag);
    switch (tag >> 3) {
      // optional .apollo.drivers.canbus.CANCardParameter.CANCardBrand brand = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 8)) {
          uint64_t val = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
          if (PROTOBUF_PREDICT_TRUE(::apollo::drivers::canbus::CANCardParameter_CANCardBrand_IsValid(val))) {
            _internal_set_brand(static_cast<::apollo::drivers::canbus::CANCardParameter_CANCardBrand>(val));
          } else {
            ::PROTOBUF_NAMESPACE_ID::internal::WriteVarint(1, val, mutable_unknown_fields());
          }
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.drivers.canbus.CANCardParameter.CANCardType type = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 16)) {
          uint64_t val = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
          if (PROTOBUF_PREDICT_TRUE(::apollo::drivers::canbus::CANCardParameter_CANCardType_IsValid(val))) {
            _internal_set_type(static_cast<::apollo::drivers::canbus::CANCardParameter_CANCardType>(val));
          } else {
            ::PROTOBUF_NAMESPACE_ID::internal::WriteVarint(2, val, mutable_unknown_fields());
          }
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.drivers.canbus.CANCardParameter.CANChannelId channel_id = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 24)) {
          uint64_t val = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
          if (PROTOBUF_PREDICT_TRUE(::apollo::drivers::canbus::CANCardParameter_CANChannelId_IsValid(val))) {
            _internal_set_channel_id(static_cast<::apollo::drivers::canbus::CANCardParameter_CANChannelId>(val));
          } else {
            ::PROTOBUF_NAMESPACE_ID::internal::WriteVarint(3, val, mutable_unknown_fields());
          }
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.drivers.canbus.CANCardParameter.CANInterface interface = 4;
      case 4:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 32)) {
          uint64_t val = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
          if (PROTOBUF_PREDICT_TRUE(::apollo::drivers::canbus::CANCardParameter_CANInterface_IsValid(val))) {
            _internal_set_interface(static_cast<::apollo::drivers::canbus::CANCardParameter_CANInterface>(val));
          } else {
            ::PROTOBUF_NAMESPACE_ID::internal::WriteVarint(4, val, mutable_unknown_fields());
          }
        } else
          goto handle_unusual;
        continue;
      // optional uint32 num_ports = 5 [default = 4];
      case 5:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 40)) {
          _Internal::set_has_num_ports(&has_bits);
          _impl_.num_ports_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint32(&ptr);
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

uint8_t* CANCardParameter::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:apollo.drivers.canbus.CANCardParameter)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  // optional .apollo.drivers.canbus.CANCardParameter.CANCardBrand brand = 1;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteEnumToArray(
      1, this->_internal_brand(), target);
  }

  // optional .apollo.drivers.canbus.CANCardParameter.CANCardType type = 2;
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteEnumToArray(
      2, this->_internal_type(), target);
  }

  // optional .apollo.drivers.canbus.CANCardParameter.CANChannelId channel_id = 3;
  if (cached_has_bits & 0x00000004u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteEnumToArray(
      3, this->_internal_channel_id(), target);
  }

  // optional .apollo.drivers.canbus.CANCardParameter.CANInterface interface = 4;
  if (cached_has_bits & 0x00000008u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteEnumToArray(
      4, this->_internal_interface(), target);
  }

  // optional uint32 num_ports = 5 [default = 4];
  if (cached_has_bits & 0x00000010u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteUInt32ToArray(5, this->_internal_num_ports(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:apollo.drivers.canbus.CANCardParameter)
  return target;
}

size_t CANCardParameter::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:apollo.drivers.canbus.CANCardParameter)
  size_t total_size = 0;

  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    // optional .apollo.drivers.canbus.CANCardParameter.CANCardBrand brand = 1;
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 +
        ::_pbi::WireFormatLite::EnumSize(this->_internal_brand());
    }

    // optional .apollo.drivers.canbus.CANCardParameter.CANCardType type = 2;
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 +
        ::_pbi::WireFormatLite::EnumSize(this->_internal_type());
    }

    // optional .apollo.drivers.canbus.CANCardParameter.CANChannelId channel_id = 3;
    if (cached_has_bits & 0x00000004u) {
      total_size += 1 +
        ::_pbi::WireFormatLite::EnumSize(this->_internal_channel_id());
    }

    // optional .apollo.drivers.canbus.CANCardParameter.CANInterface interface = 4;
    if (cached_has_bits & 0x00000008u) {
      total_size += 1 +
        ::_pbi::WireFormatLite::EnumSize(this->_internal_interface());
    }

    // optional uint32 num_ports = 5 [default = 4];
    if (cached_has_bits & 0x00000010u) {
      total_size += ::_pbi::WireFormatLite::UInt32SizePlusOne(this->_internal_num_ports());
    }

  }
  return MaybeComputeUnknownFieldsSize(total_size, &_impl_._cached_size_);
}

const ::PROTOBUF_NAMESPACE_ID::Message::ClassData CANCardParameter::_class_data_ = {
    ::PROTOBUF_NAMESPACE_ID::Message::CopyWithSizeCheck,
    CANCardParameter::MergeImpl
};
const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*CANCardParameter::GetClassData() const { return &_class_data_; }

void CANCardParameter::MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to,
                      const ::PROTOBUF_NAMESPACE_ID::Message& from) {
  static_cast<CANCardParameter *>(to)->MergeFrom(
      static_cast<const CANCardParameter &>(from));
}


void CANCardParameter::MergeFrom(const CANCardParameter& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:apollo.drivers.canbus.CANCardParameter)
  GOOGLE_DCHECK_NE(&from, this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._impl_._has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    if (cached_has_bits & 0x00000001u) {
      _impl_.brand_ = from._impl_.brand_;
    }
    if (cached_has_bits & 0x00000002u) {
      _impl_.type_ = from._impl_.type_;
    }
    if (cached_has_bits & 0x00000004u) {
      _impl_.channel_id_ = from._impl_.channel_id_;
    }
    if (cached_has_bits & 0x00000008u) {
      _impl_.interface_ = from._impl_.interface_;
    }
    if (cached_has_bits & 0x00000010u) {
      _impl_.num_ports_ = from._impl_.num_ports_;
    }
    _impl_._has_bits_[0] |= cached_has_bits;
  }
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void CANCardParameter::CopyFrom(const CANCardParameter& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:apollo.drivers.canbus.CANCardParameter)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool CANCardParameter::IsInitialized() const {
  return true;
}

void CANCardParameter::InternalSwap(CANCardParameter* other) {
  using std::swap;
  _internal_metadata_.InternalSwap(&other->_internal_metadata_);
  swap(_impl_._has_bits_[0], other->_impl_._has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(CANCardParameter, _impl_.interface_)
      + sizeof(CANCardParameter::_impl_.interface_)
      - PROTOBUF_FIELD_OFFSET(CANCardParameter, _impl_.brand_)>(
          reinterpret_cast<char*>(&_impl_.brand_),
          reinterpret_cast<char*>(&other->_impl_.brand_));
  swap(_impl_.num_ports_, other->_impl_.num_ports_);
}

::PROTOBUF_NAMESPACE_ID::Metadata CANCardParameter::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto_getter, &descriptor_table_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto_once,
      file_level_metadata_modules_2fdrivers_2fcanbus_2fproto_2fcan_5fcard_5fparameter_2eproto[0]);
}

// @@protoc_insertion_point(namespace_scope)
}  // namespace canbus
}  // namespace drivers
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::apollo::drivers::canbus::CANCardParameter*
Arena::CreateMaybeMessage< ::apollo::drivers::canbus::CANCardParameter >(Arena* arena) {
  return Arena::CreateMessageInternal< ::apollo::drivers::canbus::CANCardParameter >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
