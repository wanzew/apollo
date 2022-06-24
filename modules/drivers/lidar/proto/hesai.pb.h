// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/drivers/lidar/proto/hesai.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_modules_2fdrivers_2flidar_2fproto_2fhesai_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_modules_2fdrivers_2flidar_2fproto_2fhesai_2eproto

#include <limits>
#include <string>

#include <google/protobuf/port_def.inc>
#if PROTOBUF_VERSION < 3020000
#error This file was generated by a newer version of protoc which is
#error incompatible with your Protocol Buffer headers. Please update
#error your headers.
#endif
#if 3020000 < PROTOBUF_MIN_PROTOC_VERSION
#error This file was generated by an older version of protoc which is
#error incompatible with your Protocol Buffer headers. Please
#error regenerate this file with a newer version of protoc.
#endif

#include <google/protobuf/port_undef.inc>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/arena.h>
#include <google/protobuf/arenastring.h>
#include <google/protobuf/generated_message_util.h>
#include <google/protobuf/metadata_lite.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/message.h>
#include <google/protobuf/repeated_field.h>  // IWYU pragma: export
#include <google/protobuf/extension_set.h>  // IWYU pragma: export
#include <google/protobuf/generated_enum_reflection.h>
#include <google/protobuf/unknown_field_set.h>
#include "modules/common/proto/header.pb.h"
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_modules_2fdrivers_2flidar_2fproto_2fhesai_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_modules_2fdrivers_2flidar_2fproto_2fhesai_2eproto {
  static const uint32_t offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_modules_2fdrivers_2flidar_2fproto_2fhesai_2eproto;
namespace apollo {
namespace drivers {
namespace hesai {
class HesaiScan;
struct HesaiScanDefaultTypeInternal;
extern HesaiScanDefaultTypeInternal _HesaiScan_default_instance_;
class HesaiScanPacket;
struct HesaiScanPacketDefaultTypeInternal;
extern HesaiScanPacketDefaultTypeInternal _HesaiScanPacket_default_instance_;
}  // namespace hesai
}  // namespace drivers
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> ::apollo::drivers::hesai::HesaiScan* Arena::CreateMaybeMessage<::apollo::drivers::hesai::HesaiScan>(Arena*);
template<> ::apollo::drivers::hesai::HesaiScanPacket* Arena::CreateMaybeMessage<::apollo::drivers::hesai::HesaiScanPacket>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace apollo {
namespace drivers {
namespace hesai {

enum Model : int {
  UNKNOWN = 0,
  HESAI40P = 1,
  HESAI64 = 2
};
bool Model_IsValid(int value);
constexpr Model Model_MIN = UNKNOWN;
constexpr Model Model_MAX = HESAI64;
constexpr int Model_ARRAYSIZE = Model_MAX + 1;

const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* Model_descriptor();
template<typename T>
inline const std::string& Model_Name(T enum_t_value) {
  static_assert(::std::is_same<T, Model>::value ||
    ::std::is_integral<T>::value,
    "Incorrect type passed to function Model_Name.");
  return ::PROTOBUF_NAMESPACE_ID::internal::NameOfEnum(
    Model_descriptor(), enum_t_value);
}
inline bool Model_Parse(
    ::PROTOBUF_NAMESPACE_ID::ConstStringParam name, Model* value) {
  return ::PROTOBUF_NAMESPACE_ID::internal::ParseNamedEnum<Model>(
    Model_descriptor(), name, value);
}
// ===================================================================

class HesaiScanPacket final :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:apollo.drivers.hesai.HesaiScanPacket) */ {
 public:
  inline HesaiScanPacket() : HesaiScanPacket(nullptr) {}
  ~HesaiScanPacket() override;
  explicit PROTOBUF_CONSTEXPR HesaiScanPacket(::PROTOBUF_NAMESPACE_ID::internal::ConstantInitialized);

  HesaiScanPacket(const HesaiScanPacket& from);
  HesaiScanPacket(HesaiScanPacket&& from) noexcept
    : HesaiScanPacket() {
    *this = ::std::move(from);
  }

  inline HesaiScanPacket& operator=(const HesaiScanPacket& from) {
    CopyFrom(from);
    return *this;
  }
  inline HesaiScanPacket& operator=(HesaiScanPacket&& from) noexcept {
    if (this == &from) return *this;
    if (GetOwningArena() == from.GetOwningArena()
  #ifdef PROTOBUF_FORCE_COPY_IN_MOVE
        && GetOwningArena() != nullptr
  #endif  // !PROTOBUF_FORCE_COPY_IN_MOVE
    ) {
      InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }

  inline const ::PROTOBUF_NAMESPACE_ID::UnknownFieldSet& unknown_fields() const {
    return _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance);
  }
  inline ::PROTOBUF_NAMESPACE_ID::UnknownFieldSet* mutable_unknown_fields() {
    return _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
  }

  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* descriptor() {
    return GetDescriptor();
  }
  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* GetDescriptor() {
    return default_instance().GetMetadata().descriptor;
  }
  static const ::PROTOBUF_NAMESPACE_ID::Reflection* GetReflection() {
    return default_instance().GetMetadata().reflection;
  }
  static const HesaiScanPacket& default_instance() {
    return *internal_default_instance();
  }
  static inline const HesaiScanPacket* internal_default_instance() {
    return reinterpret_cast<const HesaiScanPacket*>(
               &_HesaiScanPacket_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(HesaiScanPacket& a, HesaiScanPacket& b) {
    a.Swap(&b);
  }
  inline void Swap(HesaiScanPacket* other) {
    if (other == this) return;
  #ifdef PROTOBUF_FORCE_COPY_IN_SWAP
    if (GetOwningArena() != nullptr &&
        GetOwningArena() == other->GetOwningArena()) {
   #else  // PROTOBUF_FORCE_COPY_IN_SWAP
    if (GetOwningArena() == other->GetOwningArena()) {
  #endif  // !PROTOBUF_FORCE_COPY_IN_SWAP
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(HesaiScanPacket* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetOwningArena() == other->GetOwningArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  HesaiScanPacket* New(::PROTOBUF_NAMESPACE_ID::Arena* arena = nullptr) const final {
    return CreateMaybeMessage<HesaiScanPacket>(arena);
  }
  using ::PROTOBUF_NAMESPACE_ID::Message::CopyFrom;
  void CopyFrom(const HesaiScanPacket& from);
  using ::PROTOBUF_NAMESPACE_ID::Message::MergeFrom;
  void MergeFrom(const HesaiScanPacket& from);
  private:
  static void MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to, const ::PROTOBUF_NAMESPACE_ID::Message& from);
  public:
  PROTOBUF_ATTRIBUTE_REINITIALIZES void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  const char* _InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) final;
  uint8_t* _InternalSerialize(
      uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const final;
  int GetCachedSize() const final { return _impl_._cached_size_.Get(); }

  private:
  void SharedCtor(::PROTOBUF_NAMESPACE_ID::Arena* arena, bool is_message_owned);
  void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(HesaiScanPacket* other);

  private:
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "apollo.drivers.hesai.HesaiScanPacket";
  }
  protected:
  explicit HesaiScanPacket(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                       bool is_message_owned = false);
  public:

  static const ClassData _class_data_;
  const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*GetClassData() const final;

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kDataFieldNumber = 2,
    kStampFieldNumber = 1,
  };
  // optional bytes data = 2;
  bool has_data() const;
  private:
  bool _internal_has_data() const;
  public:
  void clear_data();
  const std::string& data() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_data(ArgT0&& arg0, ArgT... args);
  std::string* mutable_data();
  PROTOBUF_NODISCARD std::string* release_data();
  void set_allocated_data(std::string* data);
  private:
  const std::string& _internal_data() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_data(const std::string& value);
  std::string* _internal_mutable_data();
  public:

  // optional uint64 stamp = 1;
  bool has_stamp() const;
  private:
  bool _internal_has_stamp() const;
  public:
  void clear_stamp();
  uint64_t stamp() const;
  void set_stamp(uint64_t value);
  private:
  uint64_t _internal_stamp() const;
  void _internal_set_stamp(uint64_t value);
  public:

  // @@protoc_insertion_point(class_scope:apollo.drivers.hesai.HesaiScanPacket)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  struct Impl_ {
    ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
    mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr data_;
    uint64_t stamp_;
  };
  union { Impl_ _impl_; };
  friend struct ::TableStruct_modules_2fdrivers_2flidar_2fproto_2fhesai_2eproto;
};
// -------------------------------------------------------------------

class HesaiScan final :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:apollo.drivers.hesai.HesaiScan) */ {
 public:
  inline HesaiScan() : HesaiScan(nullptr) {}
  ~HesaiScan() override;
  explicit PROTOBUF_CONSTEXPR HesaiScan(::PROTOBUF_NAMESPACE_ID::internal::ConstantInitialized);

  HesaiScan(const HesaiScan& from);
  HesaiScan(HesaiScan&& from) noexcept
    : HesaiScan() {
    *this = ::std::move(from);
  }

  inline HesaiScan& operator=(const HesaiScan& from) {
    CopyFrom(from);
    return *this;
  }
  inline HesaiScan& operator=(HesaiScan&& from) noexcept {
    if (this == &from) return *this;
    if (GetOwningArena() == from.GetOwningArena()
  #ifdef PROTOBUF_FORCE_COPY_IN_MOVE
        && GetOwningArena() != nullptr
  #endif  // !PROTOBUF_FORCE_COPY_IN_MOVE
    ) {
      InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }

  inline const ::PROTOBUF_NAMESPACE_ID::UnknownFieldSet& unknown_fields() const {
    return _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance);
  }
  inline ::PROTOBUF_NAMESPACE_ID::UnknownFieldSet* mutable_unknown_fields() {
    return _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
  }

  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* descriptor() {
    return GetDescriptor();
  }
  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* GetDescriptor() {
    return default_instance().GetMetadata().descriptor;
  }
  static const ::PROTOBUF_NAMESPACE_ID::Reflection* GetReflection() {
    return default_instance().GetMetadata().reflection;
  }
  static const HesaiScan& default_instance() {
    return *internal_default_instance();
  }
  static inline const HesaiScan* internal_default_instance() {
    return reinterpret_cast<const HesaiScan*>(
               &_HesaiScan_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    1;

  friend void swap(HesaiScan& a, HesaiScan& b) {
    a.Swap(&b);
  }
  inline void Swap(HesaiScan* other) {
    if (other == this) return;
  #ifdef PROTOBUF_FORCE_COPY_IN_SWAP
    if (GetOwningArena() != nullptr &&
        GetOwningArena() == other->GetOwningArena()) {
   #else  // PROTOBUF_FORCE_COPY_IN_SWAP
    if (GetOwningArena() == other->GetOwningArena()) {
  #endif  // !PROTOBUF_FORCE_COPY_IN_SWAP
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(HesaiScan* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetOwningArena() == other->GetOwningArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  HesaiScan* New(::PROTOBUF_NAMESPACE_ID::Arena* arena = nullptr) const final {
    return CreateMaybeMessage<HesaiScan>(arena);
  }
  using ::PROTOBUF_NAMESPACE_ID::Message::CopyFrom;
  void CopyFrom(const HesaiScan& from);
  using ::PROTOBUF_NAMESPACE_ID::Message::MergeFrom;
  void MergeFrom(const HesaiScan& from);
  private:
  static void MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to, const ::PROTOBUF_NAMESPACE_ID::Message& from);
  public:
  PROTOBUF_ATTRIBUTE_REINITIALIZES void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  const char* _InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) final;
  uint8_t* _InternalSerialize(
      uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const final;
  int GetCachedSize() const final { return _impl_._cached_size_.Get(); }

  private:
  void SharedCtor(::PROTOBUF_NAMESPACE_ID::Arena* arena, bool is_message_owned);
  void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(HesaiScan* other);

  private:
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "apollo.drivers.hesai.HesaiScan";
  }
  protected:
  explicit HesaiScan(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                       bool is_message_owned = false);
  public:

  static const ClassData _class_data_;
  const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*GetClassData() const final;

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kFiringPktsFieldNumber = 3,
    kHeaderFieldNumber = 1,
    kBasetimeFieldNumber = 4,
    kModelFieldNumber = 2,
  };
  // repeated .apollo.drivers.hesai.HesaiScanPacket firing_pkts = 3;
  int firing_pkts_size() const;
  private:
  int _internal_firing_pkts_size() const;
  public:
  void clear_firing_pkts();
  ::apollo::drivers::hesai::HesaiScanPacket* mutable_firing_pkts(int index);
  ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::apollo::drivers::hesai::HesaiScanPacket >*
      mutable_firing_pkts();
  private:
  const ::apollo::drivers::hesai::HesaiScanPacket& _internal_firing_pkts(int index) const;
  ::apollo::drivers::hesai::HesaiScanPacket* _internal_add_firing_pkts();
  public:
  const ::apollo::drivers::hesai::HesaiScanPacket& firing_pkts(int index) const;
  ::apollo::drivers::hesai::HesaiScanPacket* add_firing_pkts();
  const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::apollo::drivers::hesai::HesaiScanPacket >&
      firing_pkts() const;

  // optional .apollo.common.Header header = 1;
  bool has_header() const;
  private:
  bool _internal_has_header() const;
  public:
  void clear_header();
  const ::apollo::common::Header& header() const;
  PROTOBUF_NODISCARD ::apollo::common::Header* release_header();
  ::apollo::common::Header* mutable_header();
  void set_allocated_header(::apollo::common::Header* header);
  private:
  const ::apollo::common::Header& _internal_header() const;
  ::apollo::common::Header* _internal_mutable_header();
  public:
  void unsafe_arena_set_allocated_header(
      ::apollo::common::Header* header);
  ::apollo::common::Header* unsafe_arena_release_header();

  // optional uint64 basetime = 4 [default = 0];
  bool has_basetime() const;
  private:
  bool _internal_has_basetime() const;
  public:
  void clear_basetime();
  uint64_t basetime() const;
  void set_basetime(uint64_t value);
  private:
  uint64_t _internal_basetime() const;
  void _internal_set_basetime(uint64_t value);
  public:

  // optional .apollo.drivers.hesai.Model model = 2;
  bool has_model() const;
  private:
  bool _internal_has_model() const;
  public:
  void clear_model();
  ::apollo::drivers::hesai::Model model() const;
  void set_model(::apollo::drivers::hesai::Model value);
  private:
  ::apollo::drivers::hesai::Model _internal_model() const;
  void _internal_set_model(::apollo::drivers::hesai::Model value);
  public:

  // @@protoc_insertion_point(class_scope:apollo.drivers.hesai.HesaiScan)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  struct Impl_ {
    ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
    mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
    ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::apollo::drivers::hesai::HesaiScanPacket > firing_pkts_;
    ::apollo::common::Header* header_;
    uint64_t basetime_;
    int model_;
  };
  union { Impl_ _impl_; };
  friend struct ::TableStruct_modules_2fdrivers_2flidar_2fproto_2fhesai_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// HesaiScanPacket

// optional uint64 stamp = 1;
inline bool HesaiScanPacket::_internal_has_stamp() const {
  bool value = (_impl_._has_bits_[0] & 0x00000002u) != 0;
  return value;
}
inline bool HesaiScanPacket::has_stamp() const {
  return _internal_has_stamp();
}
inline void HesaiScanPacket::clear_stamp() {
  _impl_.stamp_ = uint64_t{0u};
  _impl_._has_bits_[0] &= ~0x00000002u;
}
inline uint64_t HesaiScanPacket::_internal_stamp() const {
  return _impl_.stamp_;
}
inline uint64_t HesaiScanPacket::stamp() const {
  // @@protoc_insertion_point(field_get:apollo.drivers.hesai.HesaiScanPacket.stamp)
  return _internal_stamp();
}
inline void HesaiScanPacket::_internal_set_stamp(uint64_t value) {
  _impl_._has_bits_[0] |= 0x00000002u;
  _impl_.stamp_ = value;
}
inline void HesaiScanPacket::set_stamp(uint64_t value) {
  _internal_set_stamp(value);
  // @@protoc_insertion_point(field_set:apollo.drivers.hesai.HesaiScanPacket.stamp)
}

// optional bytes data = 2;
inline bool HesaiScanPacket::_internal_has_data() const {
  bool value = (_impl_._has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool HesaiScanPacket::has_data() const {
  return _internal_has_data();
}
inline void HesaiScanPacket::clear_data() {
  _impl_.data_.ClearToEmpty();
  _impl_._has_bits_[0] &= ~0x00000001u;
}
inline const std::string& HesaiScanPacket::data() const {
  // @@protoc_insertion_point(field_get:apollo.drivers.hesai.HesaiScanPacket.data)
  return _internal_data();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void HesaiScanPacket::set_data(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000001u;
 _impl_.data_.SetBytes(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.drivers.hesai.HesaiScanPacket.data)
}
inline std::string* HesaiScanPacket::mutable_data() {
  std::string* _s = _internal_mutable_data();
  // @@protoc_insertion_point(field_mutable:apollo.drivers.hesai.HesaiScanPacket.data)
  return _s;
}
inline const std::string& HesaiScanPacket::_internal_data() const {
  return _impl_.data_.Get();
}
inline void HesaiScanPacket::_internal_set_data(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000001u;
  _impl_.data_.Set(value, GetArenaForAllocation());
}
inline std::string* HesaiScanPacket::_internal_mutable_data() {
  _impl_._has_bits_[0] |= 0x00000001u;
  return _impl_.data_.Mutable(GetArenaForAllocation());
}
inline std::string* HesaiScanPacket::release_data() {
  // @@protoc_insertion_point(field_release:apollo.drivers.hesai.HesaiScanPacket.data)
  if (!_internal_has_data()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000001u;
  auto* p = _impl_.data_.Release();
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.data_.IsDefault()) {
    _impl_.data_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  return p;
}
inline void HesaiScanPacket::set_allocated_data(std::string* data) {
  if (data != nullptr) {
    _impl_._has_bits_[0] |= 0x00000001u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000001u;
  }
  _impl_.data_.SetAllocated(data, GetArenaForAllocation());
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.data_.IsDefault()) {
    _impl_.data_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  // @@protoc_insertion_point(field_set_allocated:apollo.drivers.hesai.HesaiScanPacket.data)
}

// -------------------------------------------------------------------

// HesaiScan

// optional .apollo.common.Header header = 1;
inline bool HesaiScan::_internal_has_header() const {
  bool value = (_impl_._has_bits_[0] & 0x00000001u) != 0;
  PROTOBUF_ASSUME(!value || _impl_.header_ != nullptr);
  return value;
}
inline bool HesaiScan::has_header() const {
  return _internal_has_header();
}
inline const ::apollo::common::Header& HesaiScan::_internal_header() const {
  const ::apollo::common::Header* p = _impl_.header_;
  return p != nullptr ? *p : reinterpret_cast<const ::apollo::common::Header&>(
      ::apollo::common::_Header_default_instance_);
}
inline const ::apollo::common::Header& HesaiScan::header() const {
  // @@protoc_insertion_point(field_get:apollo.drivers.hesai.HesaiScan.header)
  return _internal_header();
}
inline void HesaiScan::unsafe_arena_set_allocated_header(
    ::apollo::common::Header* header) {
  if (GetArenaForAllocation() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(_impl_.header_);
  }
  _impl_.header_ = header;
  if (header) {
    _impl_._has_bits_[0] |= 0x00000001u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000001u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:apollo.drivers.hesai.HesaiScan.header)
}
inline ::apollo::common::Header* HesaiScan::release_header() {
  _impl_._has_bits_[0] &= ~0x00000001u;
  ::apollo::common::Header* temp = _impl_.header_;
  _impl_.header_ = nullptr;
#ifdef PROTOBUF_FORCE_COPY_IN_RELEASE
  auto* old =  reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(temp);
  temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  if (GetArenaForAllocation() == nullptr) { delete old; }
#else  // PROTOBUF_FORCE_COPY_IN_RELEASE
  if (GetArenaForAllocation() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
#endif  // !PROTOBUF_FORCE_COPY_IN_RELEASE
  return temp;
}
inline ::apollo::common::Header* HesaiScan::unsafe_arena_release_header() {
  // @@protoc_insertion_point(field_release:apollo.drivers.hesai.HesaiScan.header)
  _impl_._has_bits_[0] &= ~0x00000001u;
  ::apollo::common::Header* temp = _impl_.header_;
  _impl_.header_ = nullptr;
  return temp;
}
inline ::apollo::common::Header* HesaiScan::_internal_mutable_header() {
  _impl_._has_bits_[0] |= 0x00000001u;
  if (_impl_.header_ == nullptr) {
    auto* p = CreateMaybeMessage<::apollo::common::Header>(GetArenaForAllocation());
    _impl_.header_ = p;
  }
  return _impl_.header_;
}
inline ::apollo::common::Header* HesaiScan::mutable_header() {
  ::apollo::common::Header* _msg = _internal_mutable_header();
  // @@protoc_insertion_point(field_mutable:apollo.drivers.hesai.HesaiScan.header)
  return _msg;
}
inline void HesaiScan::set_allocated_header(::apollo::common::Header* header) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArenaForAllocation();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(_impl_.header_);
  }
  if (header) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
        ::PROTOBUF_NAMESPACE_ID::Arena::InternalGetOwningArena(
                reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(header));
    if (message_arena != submessage_arena) {
      header = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, header, submessage_arena);
    }
    _impl_._has_bits_[0] |= 0x00000001u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000001u;
  }
  _impl_.header_ = header;
  // @@protoc_insertion_point(field_set_allocated:apollo.drivers.hesai.HesaiScan.header)
}

// optional .apollo.drivers.hesai.Model model = 2;
inline bool HesaiScan::_internal_has_model() const {
  bool value = (_impl_._has_bits_[0] & 0x00000004u) != 0;
  return value;
}
inline bool HesaiScan::has_model() const {
  return _internal_has_model();
}
inline void HesaiScan::clear_model() {
  _impl_.model_ = 0;
  _impl_._has_bits_[0] &= ~0x00000004u;
}
inline ::apollo::drivers::hesai::Model HesaiScan::_internal_model() const {
  return static_cast< ::apollo::drivers::hesai::Model >(_impl_.model_);
}
inline ::apollo::drivers::hesai::Model HesaiScan::model() const {
  // @@protoc_insertion_point(field_get:apollo.drivers.hesai.HesaiScan.model)
  return _internal_model();
}
inline void HesaiScan::_internal_set_model(::apollo::drivers::hesai::Model value) {
  assert(::apollo::drivers::hesai::Model_IsValid(value));
  _impl_._has_bits_[0] |= 0x00000004u;
  _impl_.model_ = value;
}
inline void HesaiScan::set_model(::apollo::drivers::hesai::Model value) {
  _internal_set_model(value);
  // @@protoc_insertion_point(field_set:apollo.drivers.hesai.HesaiScan.model)
}

// repeated .apollo.drivers.hesai.HesaiScanPacket firing_pkts = 3;
inline int HesaiScan::_internal_firing_pkts_size() const {
  return _impl_.firing_pkts_.size();
}
inline int HesaiScan::firing_pkts_size() const {
  return _internal_firing_pkts_size();
}
inline void HesaiScan::clear_firing_pkts() {
  _impl_.firing_pkts_.Clear();
}
inline ::apollo::drivers::hesai::HesaiScanPacket* HesaiScan::mutable_firing_pkts(int index) {
  // @@protoc_insertion_point(field_mutable:apollo.drivers.hesai.HesaiScan.firing_pkts)
  return _impl_.firing_pkts_.Mutable(index);
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::apollo::drivers::hesai::HesaiScanPacket >*
HesaiScan::mutable_firing_pkts() {
  // @@protoc_insertion_point(field_mutable_list:apollo.drivers.hesai.HesaiScan.firing_pkts)
  return &_impl_.firing_pkts_;
}
inline const ::apollo::drivers::hesai::HesaiScanPacket& HesaiScan::_internal_firing_pkts(int index) const {
  return _impl_.firing_pkts_.Get(index);
}
inline const ::apollo::drivers::hesai::HesaiScanPacket& HesaiScan::firing_pkts(int index) const {
  // @@protoc_insertion_point(field_get:apollo.drivers.hesai.HesaiScan.firing_pkts)
  return _internal_firing_pkts(index);
}
inline ::apollo::drivers::hesai::HesaiScanPacket* HesaiScan::_internal_add_firing_pkts() {
  return _impl_.firing_pkts_.Add();
}
inline ::apollo::drivers::hesai::HesaiScanPacket* HesaiScan::add_firing_pkts() {
  ::apollo::drivers::hesai::HesaiScanPacket* _add = _internal_add_firing_pkts();
  // @@protoc_insertion_point(field_add:apollo.drivers.hesai.HesaiScan.firing_pkts)
  return _add;
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::apollo::drivers::hesai::HesaiScanPacket >&
HesaiScan::firing_pkts() const {
  // @@protoc_insertion_point(field_list:apollo.drivers.hesai.HesaiScan.firing_pkts)
  return _impl_.firing_pkts_;
}

// optional uint64 basetime = 4 [default = 0];
inline bool HesaiScan::_internal_has_basetime() const {
  bool value = (_impl_._has_bits_[0] & 0x00000002u) != 0;
  return value;
}
inline bool HesaiScan::has_basetime() const {
  return _internal_has_basetime();
}
inline void HesaiScan::clear_basetime() {
  _impl_.basetime_ = uint64_t{0u};
  _impl_._has_bits_[0] &= ~0x00000002u;
}
inline uint64_t HesaiScan::_internal_basetime() const {
  return _impl_.basetime_;
}
inline uint64_t HesaiScan::basetime() const {
  // @@protoc_insertion_point(field_get:apollo.drivers.hesai.HesaiScan.basetime)
  return _internal_basetime();
}
inline void HesaiScan::_internal_set_basetime(uint64_t value) {
  _impl_._has_bits_[0] |= 0x00000002u;
  _impl_.basetime_ = value;
}
inline void HesaiScan::set_basetime(uint64_t value) {
  _internal_set_basetime(value);
  // @@protoc_insertion_point(field_set:apollo.drivers.hesai.HesaiScan.basetime)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__
// -------------------------------------------------------------------


// @@protoc_insertion_point(namespace_scope)

}  // namespace hesai
}  // namespace drivers
}  // namespace apollo

PROTOBUF_NAMESPACE_OPEN

template <> struct is_proto_enum< ::apollo::drivers::hesai::Model> : ::std::true_type {};
template <>
inline const EnumDescriptor* GetEnumDescriptor< ::apollo::drivers::hesai::Model>() {
  return ::apollo::drivers::hesai::Model_descriptor();
}

PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_modules_2fdrivers_2flidar_2fproto_2fhesai_2eproto
