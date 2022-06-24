// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/common/util/testdata/simple.proto

#include "modules/common/util/testdata/simple.pb.h"

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
namespace common {
namespace util {
namespace test {
PROTOBUF_CONSTEXPR SimpleMessage::SimpleMessage(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_._has_bits_)*/{}
  , /*decltype(_impl_._cached_size_)*/{}
  , /*decltype(_impl_.text_)*/{&::_pbi::fixed_address_empty_string, ::_pbi::ConstantInitialized{}}
  , /*decltype(_impl_.header_)*/nullptr
  , /*decltype(_impl_.integer_)*/0} {}
struct SimpleMessageDefaultTypeInternal {
  PROTOBUF_CONSTEXPR SimpleMessageDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~SimpleMessageDefaultTypeInternal() {}
  union {
    SimpleMessage _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 SimpleMessageDefaultTypeInternal _SimpleMessage_default_instance_;
PROTOBUF_CONSTEXPR SimpleRepeatedMessage::SimpleRepeatedMessage(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_.message_)*/{}
  , /*decltype(_impl_._cached_size_)*/{}} {}
struct SimpleRepeatedMessageDefaultTypeInternal {
  PROTOBUF_CONSTEXPR SimpleRepeatedMessageDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~SimpleRepeatedMessageDefaultTypeInternal() {}
  union {
    SimpleRepeatedMessage _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 SimpleRepeatedMessageDefaultTypeInternal _SimpleRepeatedMessage_default_instance_;
}  // namespace test
}  // namespace util
}  // namespace common
}  // namespace apollo
static ::_pb::Metadata file_level_metadata_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto[2];
static constexpr ::_pb::EnumDescriptor const** file_level_enum_descriptors_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto = nullptr;
static constexpr ::_pb::ServiceDescriptor const** file_level_service_descriptors_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto = nullptr;

const uint32_t TableStruct_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::apollo::common::util::test::SimpleMessage, _impl_._has_bits_),
  PROTOBUF_FIELD_OFFSET(::apollo::common::util::test::SimpleMessage, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::apollo::common::util::test::SimpleMessage, _impl_.integer_),
  PROTOBUF_FIELD_OFFSET(::apollo::common::util::test::SimpleMessage, _impl_.text_),
  PROTOBUF_FIELD_OFFSET(::apollo::common::util::test::SimpleMessage, _impl_.header_),
  2,
  0,
  1,
  ~0u,  // no _has_bits_
  PROTOBUF_FIELD_OFFSET(::apollo::common::util::test::SimpleRepeatedMessage, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::apollo::common::util::test::SimpleRepeatedMessage, _impl_.message_),
};
static const ::_pbi::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 9, -1, sizeof(::apollo::common::util::test::SimpleMessage)},
  { 12, -1, -1, sizeof(::apollo::common::util::test::SimpleRepeatedMessage)},
};

static const ::_pb::Message* const file_default_instances[] = {
  &::apollo::common::util::test::_SimpleMessage_default_instance_._instance,
  &::apollo::common::util::test::_SimpleRepeatedMessage_default_instance_._instance,
};

const char descriptor_table_protodef_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n)modules/common/util/testdata/simple.pr"
  "oto\022\027apollo.common.util.test\032!modules/co"
  "mmon/proto/header.proto\"U\n\rSimpleMessage"
  "\022\017\n\007integer\030\001 \001(\005\022\014\n\004text\030\002 \001(\t\022%\n\006heade"
  "r\030\003 \001(\0132\025.apollo.common.Header\"P\n\025Simple"
  "RepeatedMessage\0227\n\007message\030\001 \003(\0132&.apoll"
  "o.common.util.test.SimpleMessage"
  ;
static const ::_pbi::DescriptorTable* const descriptor_table_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto_deps[1] = {
  &::descriptor_table_modules_2fcommon_2fproto_2fheader_2eproto,
};
static ::_pbi::once_flag descriptor_table_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto_once;
const ::_pbi::DescriptorTable descriptor_table_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto = {
    false, false, 272, descriptor_table_protodef_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto,
    "modules/common/util/testdata/simple.proto",
    &descriptor_table_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto_once, descriptor_table_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto_deps, 1, 2,
    schemas, file_default_instances, TableStruct_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto::offsets,
    file_level_metadata_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto, file_level_enum_descriptors_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto,
    file_level_service_descriptors_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto,
};
PROTOBUF_ATTRIBUTE_WEAK const ::_pbi::DescriptorTable* descriptor_table_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto_getter() {
  return &descriptor_table_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto;
}

// Force running AddDescriptors() at dynamic initialization time.
PROTOBUF_ATTRIBUTE_INIT_PRIORITY2 static ::_pbi::AddDescriptorsRunner dynamic_init_dummy_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto(&descriptor_table_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto);
namespace apollo {
namespace common {
namespace util {
namespace test {

// ===================================================================

class SimpleMessage::_Internal {
 public:
  using HasBits = decltype(std::declval<SimpleMessage>()._impl_._has_bits_);
  static void set_has_integer(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static void set_has_text(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static const ::apollo::common::Header& header(const SimpleMessage* msg);
  static void set_has_header(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
};

const ::apollo::common::Header&
SimpleMessage::_Internal::header(const SimpleMessage* msg) {
  return *msg->_impl_.header_;
}
void SimpleMessage::clear_header() {
  if (_impl_.header_ != nullptr) _impl_.header_->Clear();
  _impl_._has_bits_[0] &= ~0x00000002u;
}
SimpleMessage::SimpleMessage(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:apollo.common.util.test.SimpleMessage)
}
SimpleMessage::SimpleMessage(const SimpleMessage& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){from._impl_._has_bits_}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.text_){}
    , decltype(_impl_.header_){nullptr}
    , decltype(_impl_.integer_){}};

  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  _impl_.text_.InitDefault();
  #ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
    _impl_.text_.Set("", GetArenaForAllocation());
  #endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (from._internal_has_text()) {
    _impl_.text_.Set(from._internal_text(), 
      GetArenaForAllocation());
  }
  if (from._internal_has_header()) {
    _impl_.header_ = new ::apollo::common::Header(*from._impl_.header_);
  }
  _impl_.integer_ = from._impl_.integer_;
  // @@protoc_insertion_point(copy_constructor:apollo.common.util.test.SimpleMessage)
}

inline void SimpleMessage::SharedCtor(
    ::_pb::Arena* arena, bool is_message_owned) {
  (void)arena;
  (void)is_message_owned;
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.text_){}
    , decltype(_impl_.header_){nullptr}
    , decltype(_impl_.integer_){0}
  };
  _impl_.text_.InitDefault();
  #ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
    _impl_.text_.Set("", GetArenaForAllocation());
  #endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
}

SimpleMessage::~SimpleMessage() {
  // @@protoc_insertion_point(destructor:apollo.common.util.test.SimpleMessage)
  if (auto *arena = _internal_metadata_.DeleteReturnArena<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>()) {
  (void)arena;
    return;
  }
  SharedDtor();
}

inline void SimpleMessage::SharedDtor() {
  GOOGLE_DCHECK(GetArenaForAllocation() == nullptr);
  _impl_.text_.Destroy();
  if (this != internal_default_instance()) delete _impl_.header_;
}

void SimpleMessage::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void SimpleMessage::Clear() {
// @@protoc_insertion_point(message_clear_start:apollo.common.util.test.SimpleMessage)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    if (cached_has_bits & 0x00000001u) {
      _impl_.text_.ClearNonDefaultToEmpty();
    }
    if (cached_has_bits & 0x00000002u) {
      GOOGLE_DCHECK(_impl_.header_ != nullptr);
      _impl_.header_->Clear();
    }
  }
  _impl_.integer_ = 0;
  _impl_._has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* SimpleMessage::_InternalParse(const char* ptr, ::_pbi::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  while (!ctx->Done(&ptr)) {
    uint32_t tag;
    ptr = ::_pbi::ReadTag(ptr, &tag);
    switch (tag >> 3) {
      // optional int32 integer = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 8)) {
          _Internal::set_has_integer(&has_bits);
          _impl_.integer_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint32(&ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional string text = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 18)) {
          auto str = _internal_mutable_text();
          ptr = ::_pbi::InlineGreedyStringParser(str, ptr, ctx);
          CHK_(ptr);
          #ifndef NDEBUG
          ::_pbi::VerifyUTF8(str, "apollo.common.util.test.SimpleMessage.text");
          #endif  // !NDEBUG
        } else
          goto handle_unusual;
        continue;
      // optional .apollo.common.Header header = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 26)) {
          ptr = ctx->ParseMessage(_internal_mutable_header(), ptr);
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

uint8_t* SimpleMessage::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:apollo.common.util.test.SimpleMessage)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  // optional int32 integer = 1;
  if (cached_has_bits & 0x00000004u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteInt32ToArray(1, this->_internal_integer(), target);
  }

  // optional string text = 2;
  if (cached_has_bits & 0x00000001u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_text().data(), static_cast<int>(this->_internal_text().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "apollo.common.util.test.SimpleMessage.text");
    target = stream->WriteStringMaybeAliased(
        2, this->_internal_text(), target);
  }

  // optional .apollo.common.Header header = 3;
  if (cached_has_bits & 0x00000002u) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(3, _Internal::header(this),
        _Internal::header(this).GetCachedSize(), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:apollo.common.util.test.SimpleMessage)
  return target;
}

size_t SimpleMessage::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:apollo.common.util.test.SimpleMessage)
  size_t total_size = 0;

  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x00000007u) {
    // optional string text = 2;
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
          this->_internal_text());
    }

    // optional .apollo.common.Header header = 3;
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *_impl_.header_);
    }

    // optional int32 integer = 1;
    if (cached_has_bits & 0x00000004u) {
      total_size += ::_pbi::WireFormatLite::Int32SizePlusOne(this->_internal_integer());
    }

  }
  return MaybeComputeUnknownFieldsSize(total_size, &_impl_._cached_size_);
}

const ::PROTOBUF_NAMESPACE_ID::Message::ClassData SimpleMessage::_class_data_ = {
    ::PROTOBUF_NAMESPACE_ID::Message::CopyWithSizeCheck,
    SimpleMessage::MergeImpl
};
const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*SimpleMessage::GetClassData() const { return &_class_data_; }

void SimpleMessage::MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to,
                      const ::PROTOBUF_NAMESPACE_ID::Message& from) {
  static_cast<SimpleMessage *>(to)->MergeFrom(
      static_cast<const SimpleMessage &>(from));
}


void SimpleMessage::MergeFrom(const SimpleMessage& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:apollo.common.util.test.SimpleMessage)
  GOOGLE_DCHECK_NE(&from, this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._impl_._has_bits_[0];
  if (cached_has_bits & 0x00000007u) {
    if (cached_has_bits & 0x00000001u) {
      _internal_set_text(from._internal_text());
    }
    if (cached_has_bits & 0x00000002u) {
      _internal_mutable_header()->::apollo::common::Header::MergeFrom(from._internal_header());
    }
    if (cached_has_bits & 0x00000004u) {
      _impl_.integer_ = from._impl_.integer_;
    }
    _impl_._has_bits_[0] |= cached_has_bits;
  }
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void SimpleMessage::CopyFrom(const SimpleMessage& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:apollo.common.util.test.SimpleMessage)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool SimpleMessage::IsInitialized() const {
  return true;
}

void SimpleMessage::InternalSwap(SimpleMessage* other) {
  using std::swap;
  auto* lhs_arena = GetArenaForAllocation();
  auto* rhs_arena = other->GetArenaForAllocation();
  _internal_metadata_.InternalSwap(&other->_internal_metadata_);
  swap(_impl_._has_bits_[0], other->_impl_._has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr::InternalSwap(
      &_impl_.text_, lhs_arena,
      &other->_impl_.text_, rhs_arena
  );
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(SimpleMessage, _impl_.integer_)
      + sizeof(SimpleMessage::_impl_.integer_)
      - PROTOBUF_FIELD_OFFSET(SimpleMessage, _impl_.header_)>(
          reinterpret_cast<char*>(&_impl_.header_),
          reinterpret_cast<char*>(&other->_impl_.header_));
}

::PROTOBUF_NAMESPACE_ID::Metadata SimpleMessage::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto_getter, &descriptor_table_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto_once,
      file_level_metadata_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto[0]);
}

// ===================================================================

class SimpleRepeatedMessage::_Internal {
 public:
};

SimpleRepeatedMessage::SimpleRepeatedMessage(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:apollo.common.util.test.SimpleRepeatedMessage)
}
SimpleRepeatedMessage::SimpleRepeatedMessage(const SimpleRepeatedMessage& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  new (&_impl_) Impl_{
      decltype(_impl_.message_){from._impl_.message_}
    , /*decltype(_impl_._cached_size_)*/{}};

  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  // @@protoc_insertion_point(copy_constructor:apollo.common.util.test.SimpleRepeatedMessage)
}

inline void SimpleRepeatedMessage::SharedCtor(
    ::_pb::Arena* arena, bool is_message_owned) {
  (void)arena;
  (void)is_message_owned;
  new (&_impl_) Impl_{
      decltype(_impl_.message_){arena}
    , /*decltype(_impl_._cached_size_)*/{}
  };
}

SimpleRepeatedMessage::~SimpleRepeatedMessage() {
  // @@protoc_insertion_point(destructor:apollo.common.util.test.SimpleRepeatedMessage)
  if (auto *arena = _internal_metadata_.DeleteReturnArena<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>()) {
  (void)arena;
    return;
  }
  SharedDtor();
}

inline void SimpleRepeatedMessage::SharedDtor() {
  GOOGLE_DCHECK(GetArenaForAllocation() == nullptr);
  _impl_.message_.~RepeatedPtrField();
}

void SimpleRepeatedMessage::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void SimpleRepeatedMessage::Clear() {
// @@protoc_insertion_point(message_clear_start:apollo.common.util.test.SimpleRepeatedMessage)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  _impl_.message_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* SimpleRepeatedMessage::_InternalParse(const char* ptr, ::_pbi::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  while (!ctx->Done(&ptr)) {
    uint32_t tag;
    ptr = ::_pbi::ReadTag(ptr, &tag);
    switch (tag >> 3) {
      // repeated .apollo.common.util.test.SimpleMessage message = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 10)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_message(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<10>(ptr));
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
  return ptr;
failure:
  ptr = nullptr;
  goto message_done;
#undef CHK_
}

uint8_t* SimpleRepeatedMessage::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:apollo.common.util.test.SimpleRepeatedMessage)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated .apollo.common.util.test.SimpleMessage message = 1;
  for (unsigned i = 0,
      n = static_cast<unsigned>(this->_internal_message_size()); i < n; i++) {
    const auto& repfield = this->_internal_message(i);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
        InternalWriteMessage(1, repfield, repfield.GetCachedSize(), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:apollo.common.util.test.SimpleRepeatedMessage)
  return target;
}

size_t SimpleRepeatedMessage::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:apollo.common.util.test.SimpleRepeatedMessage)
  size_t total_size = 0;

  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated .apollo.common.util.test.SimpleMessage message = 1;
  total_size += 1UL * this->_internal_message_size();
  for (const auto& msg : this->_impl_.message_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  return MaybeComputeUnknownFieldsSize(total_size, &_impl_._cached_size_);
}

const ::PROTOBUF_NAMESPACE_ID::Message::ClassData SimpleRepeatedMessage::_class_data_ = {
    ::PROTOBUF_NAMESPACE_ID::Message::CopyWithSizeCheck,
    SimpleRepeatedMessage::MergeImpl
};
const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*SimpleRepeatedMessage::GetClassData() const { return &_class_data_; }

void SimpleRepeatedMessage::MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to,
                      const ::PROTOBUF_NAMESPACE_ID::Message& from) {
  static_cast<SimpleRepeatedMessage *>(to)->MergeFrom(
      static_cast<const SimpleRepeatedMessage &>(from));
}


void SimpleRepeatedMessage::MergeFrom(const SimpleRepeatedMessage& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:apollo.common.util.test.SimpleRepeatedMessage)
  GOOGLE_DCHECK_NE(&from, this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  _impl_.message_.MergeFrom(from._impl_.message_);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void SimpleRepeatedMessage::CopyFrom(const SimpleRepeatedMessage& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:apollo.common.util.test.SimpleRepeatedMessage)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool SimpleRepeatedMessage::IsInitialized() const {
  return true;
}

void SimpleRepeatedMessage::InternalSwap(SimpleRepeatedMessage* other) {
  using std::swap;
  _internal_metadata_.InternalSwap(&other->_internal_metadata_);
  _impl_.message_.InternalSwap(&other->_impl_.message_);
}

::PROTOBUF_NAMESPACE_ID::Metadata SimpleRepeatedMessage::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto_getter, &descriptor_table_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto_once,
      file_level_metadata_modules_2fcommon_2futil_2ftestdata_2fsimple_2eproto[1]);
}

// @@protoc_insertion_point(namespace_scope)
}  // namespace test
}  // namespace util
}  // namespace common
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::apollo::common::util::test::SimpleMessage*
Arena::CreateMaybeMessage< ::apollo::common::util::test::SimpleMessage >(Arena* arena) {
  return Arena::CreateMessageInternal< ::apollo::common::util::test::SimpleMessage >(arena);
}
template<> PROTOBUF_NOINLINE ::apollo::common::util::test::SimpleRepeatedMessage*
Arena::CreateMaybeMessage< ::apollo::common::util::test::SimpleRepeatedMessage >(Arena* arena) {
  return Arena::CreateMessageInternal< ::apollo::common::util::test::SimpleRepeatedMessage >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
