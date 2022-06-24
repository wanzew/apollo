// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/v2x/proto/v2x_service_obu_to_car.proto

#include "modules/v2x/proto/v2x_service_obu_to_car.pb.h"

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
namespace v2x {
PROTOBUF_CONSTEXPR StatusResponse::StatusResponse(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_._has_bits_)*/{}
  , /*decltype(_impl_._cached_size_)*/{}
  , /*decltype(_impl_.info_)*/{&::_pbi::fixed_address_empty_string, ::_pbi::ConstantInitialized{}}
  , /*decltype(_impl_.error_code_)*/int64_t{0}
  , /*decltype(_impl_.status_)*/false} {}
struct StatusResponseDefaultTypeInternal {
  PROTOBUF_CONSTEXPR StatusResponseDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~StatusResponseDefaultTypeInternal() {}
  union {
    StatusResponse _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 StatusResponseDefaultTypeInternal _StatusResponse_default_instance_;
}  // namespace v2x
}  // namespace apollo
static ::_pb::Metadata file_level_metadata_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto[1];
static constexpr ::_pb::EnumDescriptor const** file_level_enum_descriptors_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto = nullptr;
static constexpr ::_pb::ServiceDescriptor const** file_level_service_descriptors_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto = nullptr;

const uint32_t TableStruct_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::apollo::v2x::StatusResponse, _impl_._has_bits_),
  PROTOBUF_FIELD_OFFSET(::apollo::v2x::StatusResponse, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::apollo::v2x::StatusResponse, _impl_.status_),
  PROTOBUF_FIELD_OFFSET(::apollo::v2x::StatusResponse, _impl_.info_),
  PROTOBUF_FIELD_OFFSET(::apollo::v2x::StatusResponse, _impl_.error_code_),
  2,
  0,
  1,
};
static const ::_pbi::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 9, -1, sizeof(::apollo::v2x::StatusResponse)},
};

static const ::_pb::Message* const file_default_instances[] = {
  &::apollo::v2x::_StatusResponse_default_instance_._instance,
};

const char descriptor_table_protodef_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n.modules/v2x/proto/v2x_service_obu_to_c"
  "ar.proto\022\napollo.v2x\032-modules/v2x/proto/"
  "v2x_obu_traffic_light.proto\032#modules/v2x"
  "/proto/v2x_monitor.proto\032#modules/v2x/pr"
  "oto/v2x_obu_rsi.proto\032%modules/v2x/proto"
  "/v2x_obstacles.proto\"I\n\016StatusResponse\022\025"
  "\n\006status\030\001 \002(\010:\005false\022\014\n\004info\030\002 \001(\t\022\022\n\ne"
  "rror_code\030\003 \001(\0032\273\002\n\010ObuToCar\022Q\n\027SendPerc"
  "eptionObstacles\022\030.apollo.v2x.V2XObstacle"
  "s\032\032.apollo.v2x.StatusResponse\"\000\022T\n\023SendV"
  "2xTrafficLight\022\037.apollo.v2x.obu.ObuTraff"
  "icLight\032\032.apollo.v2x.StatusResponse\"\000\022B\n"
  "\nSendV2xRSI\022\026.apollo.v2x.obu.ObuRsi\032\032.ap"
  "ollo.v2x.StatusResponse\"\000\022B\n\014SendObuAlar"
  "m\022\024.apollo.v2x.ObuAlarm\032\032.apollo.v2x.Sta"
  "tusResponse\"\000"
  ;
static const ::_pbi::DescriptorTable* const descriptor_table_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto_deps[4] = {
  &::descriptor_table_modules_2fv2x_2fproto_2fv2x_5fmonitor_2eproto,
  &::descriptor_table_modules_2fv2x_2fproto_2fv2x_5fobstacles_2eproto,
  &::descriptor_table_modules_2fv2x_2fproto_2fv2x_5fobu_5frsi_2eproto,
  &::descriptor_table_modules_2fv2x_2fproto_2fv2x_5fobu_5ftraffic_5flight_2eproto,
};
static ::_pbi::once_flag descriptor_table_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto_once;
const ::_pbi::DescriptorTable descriptor_table_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto = {
    false, false, 613, descriptor_table_protodef_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto,
    "modules/v2x/proto/v2x_service_obu_to_car.proto",
    &descriptor_table_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto_once, descriptor_table_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto_deps, 4, 1,
    schemas, file_default_instances, TableStruct_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto::offsets,
    file_level_metadata_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto, file_level_enum_descriptors_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto,
    file_level_service_descriptors_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto,
};
PROTOBUF_ATTRIBUTE_WEAK const ::_pbi::DescriptorTable* descriptor_table_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto_getter() {
  return &descriptor_table_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto;
}

// Force running AddDescriptors() at dynamic initialization time.
PROTOBUF_ATTRIBUTE_INIT_PRIORITY2 static ::_pbi::AddDescriptorsRunner dynamic_init_dummy_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto(&descriptor_table_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto);
namespace apollo {
namespace v2x {

// ===================================================================

class StatusResponse::_Internal {
 public:
  using HasBits = decltype(std::declval<StatusResponse>()._impl_._has_bits_);
  static void set_has_status(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static void set_has_info(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_error_code(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000004) ^ 0x00000004) != 0;
  }
};

StatusResponse::StatusResponse(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:apollo.v2x.StatusResponse)
}
StatusResponse::StatusResponse(const StatusResponse& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){from._impl_._has_bits_}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.info_){}
    , decltype(_impl_.error_code_){}
    , decltype(_impl_.status_){}};

  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  _impl_.info_.InitDefault();
  #ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
    _impl_.info_.Set("", GetArenaForAllocation());
  #endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (from._internal_has_info()) {
    _impl_.info_.Set(from._internal_info(), 
      GetArenaForAllocation());
  }
  ::memcpy(&_impl_.error_code_, &from._impl_.error_code_,
    static_cast<size_t>(reinterpret_cast<char*>(&_impl_.status_) -
    reinterpret_cast<char*>(&_impl_.error_code_)) + sizeof(_impl_.status_));
  // @@protoc_insertion_point(copy_constructor:apollo.v2x.StatusResponse)
}

inline void StatusResponse::SharedCtor(
    ::_pb::Arena* arena, bool is_message_owned) {
  (void)arena;
  (void)is_message_owned;
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.info_){}
    , decltype(_impl_.error_code_){int64_t{0}}
    , decltype(_impl_.status_){false}
  };
  _impl_.info_.InitDefault();
  #ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
    _impl_.info_.Set("", GetArenaForAllocation());
  #endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
}

StatusResponse::~StatusResponse() {
  // @@protoc_insertion_point(destructor:apollo.v2x.StatusResponse)
  if (auto *arena = _internal_metadata_.DeleteReturnArena<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>()) {
  (void)arena;
    return;
  }
  SharedDtor();
}

inline void StatusResponse::SharedDtor() {
  GOOGLE_DCHECK(GetArenaForAllocation() == nullptr);
  _impl_.info_.Destroy();
}

void StatusResponse::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void StatusResponse::Clear() {
// @@protoc_insertion_point(message_clear_start:apollo.v2x.StatusResponse)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x00000001u) {
    _impl_.info_.ClearNonDefaultToEmpty();
  }
  if (cached_has_bits & 0x00000006u) {
    ::memset(&_impl_.error_code_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&_impl_.status_) -
        reinterpret_cast<char*>(&_impl_.error_code_)) + sizeof(_impl_.status_));
  }
  _impl_._has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* StatusResponse::_InternalParse(const char* ptr, ::_pbi::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  while (!ctx->Done(&ptr)) {
    uint32_t tag;
    ptr = ::_pbi::ReadTag(ptr, &tag);
    switch (tag >> 3) {
      // required bool status = 1 [default = false];
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 8)) {
          _Internal::set_has_status(&has_bits);
          _impl_.status_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional string info = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 18)) {
          auto str = _internal_mutable_info();
          ptr = ::_pbi::InlineGreedyStringParser(str, ptr, ctx);
          CHK_(ptr);
          #ifndef NDEBUG
          ::_pbi::VerifyUTF8(str, "apollo.v2x.StatusResponse.info");
          #endif  // !NDEBUG
        } else
          goto handle_unusual;
        continue;
      // optional int64 error_code = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 24)) {
          _Internal::set_has_error_code(&has_bits);
          _impl_.error_code_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
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

uint8_t* StatusResponse::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:apollo.v2x.StatusResponse)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  // required bool status = 1 [default = false];
  if (cached_has_bits & 0x00000004u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteBoolToArray(1, this->_internal_status(), target);
  }

  // optional string info = 2;
  if (cached_has_bits & 0x00000001u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_info().data(), static_cast<int>(this->_internal_info().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "apollo.v2x.StatusResponse.info");
    target = stream->WriteStringMaybeAliased(
        2, this->_internal_info(), target);
  }

  // optional int64 error_code = 3;
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteInt64ToArray(3, this->_internal_error_code(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:apollo.v2x.StatusResponse)
  return target;
}

size_t StatusResponse::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:apollo.v2x.StatusResponse)
  size_t total_size = 0;

  // required bool status = 1 [default = false];
  if (_internal_has_status()) {
    total_size += 1 + 1;
  }
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    // optional string info = 2;
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
          this->_internal_info());
    }

    // optional int64 error_code = 3;
    if (cached_has_bits & 0x00000002u) {
      total_size += ::_pbi::WireFormatLite::Int64SizePlusOne(this->_internal_error_code());
    }

  }
  return MaybeComputeUnknownFieldsSize(total_size, &_impl_._cached_size_);
}

const ::PROTOBUF_NAMESPACE_ID::Message::ClassData StatusResponse::_class_data_ = {
    ::PROTOBUF_NAMESPACE_ID::Message::CopyWithSizeCheck,
    StatusResponse::MergeImpl
};
const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*StatusResponse::GetClassData() const { return &_class_data_; }

void StatusResponse::MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to,
                      const ::PROTOBUF_NAMESPACE_ID::Message& from) {
  static_cast<StatusResponse *>(to)->MergeFrom(
      static_cast<const StatusResponse &>(from));
}


void StatusResponse::MergeFrom(const StatusResponse& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:apollo.v2x.StatusResponse)
  GOOGLE_DCHECK_NE(&from, this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._impl_._has_bits_[0];
  if (cached_has_bits & 0x00000007u) {
    if (cached_has_bits & 0x00000001u) {
      _internal_set_info(from._internal_info());
    }
    if (cached_has_bits & 0x00000002u) {
      _impl_.error_code_ = from._impl_.error_code_;
    }
    if (cached_has_bits & 0x00000004u) {
      _impl_.status_ = from._impl_.status_;
    }
    _impl_._has_bits_[0] |= cached_has_bits;
  }
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void StatusResponse::CopyFrom(const StatusResponse& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:apollo.v2x.StatusResponse)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool StatusResponse::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_impl_._has_bits_)) return false;
  return true;
}

void StatusResponse::InternalSwap(StatusResponse* other) {
  using std::swap;
  auto* lhs_arena = GetArenaForAllocation();
  auto* rhs_arena = other->GetArenaForAllocation();
  _internal_metadata_.InternalSwap(&other->_internal_metadata_);
  swap(_impl_._has_bits_[0], other->_impl_._has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr::InternalSwap(
      &_impl_.info_, lhs_arena,
      &other->_impl_.info_, rhs_arena
  );
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(StatusResponse, _impl_.status_)
      + sizeof(StatusResponse::_impl_.status_)
      - PROTOBUF_FIELD_OFFSET(StatusResponse, _impl_.error_code_)>(
          reinterpret_cast<char*>(&_impl_.error_code_),
          reinterpret_cast<char*>(&other->_impl_.error_code_));
}

::PROTOBUF_NAMESPACE_ID::Metadata StatusResponse::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto_getter, &descriptor_table_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto_once,
      file_level_metadata_modules_2fv2x_2fproto_2fv2x_5fservice_5fobu_5fto_5fcar_2eproto[0]);
}

// @@protoc_insertion_point(namespace_scope)
}  // namespace v2x
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::apollo::v2x::StatusResponse*
Arena::CreateMaybeMessage< ::apollo::v2x::StatusResponse >(Arena* arena) {
  return Arena::CreateMessageInternal< ::apollo::v2x::StatusResponse >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
