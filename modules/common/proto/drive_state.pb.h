// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/common/proto/drive_state.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_modules_2fcommon_2fproto_2fdrive_5fstate_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_modules_2fcommon_2fproto_2fdrive_5fstate_2eproto

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
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_modules_2fcommon_2fproto_2fdrive_5fstate_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_modules_2fcommon_2fproto_2fdrive_5fstate_2eproto {
  static const uint32_t offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_modules_2fcommon_2fproto_2fdrive_5fstate_2eproto;
namespace apollo {
namespace common {
class EngageAdvice;
struct EngageAdviceDefaultTypeInternal;
extern EngageAdviceDefaultTypeInternal _EngageAdvice_default_instance_;
}  // namespace common
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> ::apollo::common::EngageAdvice* Arena::CreateMaybeMessage<::apollo::common::EngageAdvice>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace apollo {
namespace common {

enum EngageAdvice_Advice : int {
  EngageAdvice_Advice_UNKNOWN = 0,
  EngageAdvice_Advice_DISALLOW_ENGAGE = 1,
  EngageAdvice_Advice_READY_TO_ENGAGE = 2,
  EngageAdvice_Advice_KEEP_ENGAGED = 3,
  EngageAdvice_Advice_PREPARE_DISENGAGE = 4
};
bool EngageAdvice_Advice_IsValid(int value);
constexpr EngageAdvice_Advice EngageAdvice_Advice_Advice_MIN = EngageAdvice_Advice_UNKNOWN;
constexpr EngageAdvice_Advice EngageAdvice_Advice_Advice_MAX = EngageAdvice_Advice_PREPARE_DISENGAGE;
constexpr int EngageAdvice_Advice_Advice_ARRAYSIZE = EngageAdvice_Advice_Advice_MAX + 1;

const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* EngageAdvice_Advice_descriptor();
template<typename T>
inline const std::string& EngageAdvice_Advice_Name(T enum_t_value) {
  static_assert(::std::is_same<T, EngageAdvice_Advice>::value ||
    ::std::is_integral<T>::value,
    "Incorrect type passed to function EngageAdvice_Advice_Name.");
  return ::PROTOBUF_NAMESPACE_ID::internal::NameOfEnum(
    EngageAdvice_Advice_descriptor(), enum_t_value);
}
inline bool EngageAdvice_Advice_Parse(
    ::PROTOBUF_NAMESPACE_ID::ConstStringParam name, EngageAdvice_Advice* value) {
  return ::PROTOBUF_NAMESPACE_ID::internal::ParseNamedEnum<EngageAdvice_Advice>(
    EngageAdvice_Advice_descriptor(), name, value);
}
// ===================================================================

class EngageAdvice final :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:apollo.common.EngageAdvice) */ {
 public:
  inline EngageAdvice() : EngageAdvice(nullptr) {}
  ~EngageAdvice() override;
  explicit PROTOBUF_CONSTEXPR EngageAdvice(::PROTOBUF_NAMESPACE_ID::internal::ConstantInitialized);

  EngageAdvice(const EngageAdvice& from);
  EngageAdvice(EngageAdvice&& from) noexcept
    : EngageAdvice() {
    *this = ::std::move(from);
  }

  inline EngageAdvice& operator=(const EngageAdvice& from) {
    CopyFrom(from);
    return *this;
  }
  inline EngageAdvice& operator=(EngageAdvice&& from) noexcept {
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
  static const EngageAdvice& default_instance() {
    return *internal_default_instance();
  }
  static inline const EngageAdvice* internal_default_instance() {
    return reinterpret_cast<const EngageAdvice*>(
               &_EngageAdvice_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(EngageAdvice& a, EngageAdvice& b) {
    a.Swap(&b);
  }
  inline void Swap(EngageAdvice* other) {
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
  void UnsafeArenaSwap(EngageAdvice* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetOwningArena() == other->GetOwningArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  EngageAdvice* New(::PROTOBUF_NAMESPACE_ID::Arena* arena = nullptr) const final {
    return CreateMaybeMessage<EngageAdvice>(arena);
  }
  using ::PROTOBUF_NAMESPACE_ID::Message::CopyFrom;
  void CopyFrom(const EngageAdvice& from);
  using ::PROTOBUF_NAMESPACE_ID::Message::MergeFrom;
  void MergeFrom(const EngageAdvice& from);
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
  void InternalSwap(EngageAdvice* other);

  private:
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "apollo.common.EngageAdvice";
  }
  protected:
  explicit EngageAdvice(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                       bool is_message_owned = false);
  public:

  static const ClassData _class_data_;
  const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*GetClassData() const final;

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  typedef EngageAdvice_Advice Advice;
  static constexpr Advice UNKNOWN =
    EngageAdvice_Advice_UNKNOWN;
  static constexpr Advice DISALLOW_ENGAGE =
    EngageAdvice_Advice_DISALLOW_ENGAGE;
  static constexpr Advice READY_TO_ENGAGE =
    EngageAdvice_Advice_READY_TO_ENGAGE;
  static constexpr Advice KEEP_ENGAGED =
    EngageAdvice_Advice_KEEP_ENGAGED;
  static constexpr Advice PREPARE_DISENGAGE =
    EngageAdvice_Advice_PREPARE_DISENGAGE;
  static inline bool Advice_IsValid(int value) {
    return EngageAdvice_Advice_IsValid(value);
  }
  static constexpr Advice Advice_MIN =
    EngageAdvice_Advice_Advice_MIN;
  static constexpr Advice Advice_MAX =
    EngageAdvice_Advice_Advice_MAX;
  static constexpr int Advice_ARRAYSIZE =
    EngageAdvice_Advice_Advice_ARRAYSIZE;
  static inline const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor*
  Advice_descriptor() {
    return EngageAdvice_Advice_descriptor();
  }
  template<typename T>
  static inline const std::string& Advice_Name(T enum_t_value) {
    static_assert(::std::is_same<T, Advice>::value ||
      ::std::is_integral<T>::value,
      "Incorrect type passed to function Advice_Name.");
    return EngageAdvice_Advice_Name(enum_t_value);
  }
  static inline bool Advice_Parse(::PROTOBUF_NAMESPACE_ID::ConstStringParam name,
      Advice* value) {
    return EngageAdvice_Advice_Parse(name, value);
  }

  // accessors -------------------------------------------------------

  enum : int {
    kReasonFieldNumber = 2,
    kAdviceFieldNumber = 1,
  };
  // optional string reason = 2;
  bool has_reason() const;
  private:
  bool _internal_has_reason() const;
  public:
  void clear_reason();
  const std::string& reason() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_reason(ArgT0&& arg0, ArgT... args);
  std::string* mutable_reason();
  PROTOBUF_NODISCARD std::string* release_reason();
  void set_allocated_reason(std::string* reason);
  private:
  const std::string& _internal_reason() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_reason(const std::string& value);
  std::string* _internal_mutable_reason();
  public:

  // optional .apollo.common.EngageAdvice.Advice advice = 1 [default = DISALLOW_ENGAGE];
  bool has_advice() const;
  private:
  bool _internal_has_advice() const;
  public:
  void clear_advice();
  ::apollo::common::EngageAdvice_Advice advice() const;
  void set_advice(::apollo::common::EngageAdvice_Advice value);
  private:
  ::apollo::common::EngageAdvice_Advice _internal_advice() const;
  void _internal_set_advice(::apollo::common::EngageAdvice_Advice value);
  public:

  // @@protoc_insertion_point(class_scope:apollo.common.EngageAdvice)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  struct Impl_ {
    ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
    mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr reason_;
    int advice_;
  };
  union { Impl_ _impl_; };
  friend struct ::TableStruct_modules_2fcommon_2fproto_2fdrive_5fstate_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// EngageAdvice

// optional .apollo.common.EngageAdvice.Advice advice = 1 [default = DISALLOW_ENGAGE];
inline bool EngageAdvice::_internal_has_advice() const {
  bool value = (_impl_._has_bits_[0] & 0x00000002u) != 0;
  return value;
}
inline bool EngageAdvice::has_advice() const {
  return _internal_has_advice();
}
inline void EngageAdvice::clear_advice() {
  _impl_.advice_ = 1;
  _impl_._has_bits_[0] &= ~0x00000002u;
}
inline ::apollo::common::EngageAdvice_Advice EngageAdvice::_internal_advice() const {
  return static_cast< ::apollo::common::EngageAdvice_Advice >(_impl_.advice_);
}
inline ::apollo::common::EngageAdvice_Advice EngageAdvice::advice() const {
  // @@protoc_insertion_point(field_get:apollo.common.EngageAdvice.advice)
  return _internal_advice();
}
inline void EngageAdvice::_internal_set_advice(::apollo::common::EngageAdvice_Advice value) {
  assert(::apollo::common::EngageAdvice_Advice_IsValid(value));
  _impl_._has_bits_[0] |= 0x00000002u;
  _impl_.advice_ = value;
}
inline void EngageAdvice::set_advice(::apollo::common::EngageAdvice_Advice value) {
  _internal_set_advice(value);
  // @@protoc_insertion_point(field_set:apollo.common.EngageAdvice.advice)
}

// optional string reason = 2;
inline bool EngageAdvice::_internal_has_reason() const {
  bool value = (_impl_._has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool EngageAdvice::has_reason() const {
  return _internal_has_reason();
}
inline void EngageAdvice::clear_reason() {
  _impl_.reason_.ClearToEmpty();
  _impl_._has_bits_[0] &= ~0x00000001u;
}
inline const std::string& EngageAdvice::reason() const {
  // @@protoc_insertion_point(field_get:apollo.common.EngageAdvice.reason)
  return _internal_reason();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void EngageAdvice::set_reason(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000001u;
 _impl_.reason_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.common.EngageAdvice.reason)
}
inline std::string* EngageAdvice::mutable_reason() {
  std::string* _s = _internal_mutable_reason();
  // @@protoc_insertion_point(field_mutable:apollo.common.EngageAdvice.reason)
  return _s;
}
inline const std::string& EngageAdvice::_internal_reason() const {
  return _impl_.reason_.Get();
}
inline void EngageAdvice::_internal_set_reason(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000001u;
  _impl_.reason_.Set(value, GetArenaForAllocation());
}
inline std::string* EngageAdvice::_internal_mutable_reason() {
  _impl_._has_bits_[0] |= 0x00000001u;
  return _impl_.reason_.Mutable(GetArenaForAllocation());
}
inline std::string* EngageAdvice::release_reason() {
  // @@protoc_insertion_point(field_release:apollo.common.EngageAdvice.reason)
  if (!_internal_has_reason()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000001u;
  auto* p = _impl_.reason_.Release();
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.reason_.IsDefault()) {
    _impl_.reason_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  return p;
}
inline void EngageAdvice::set_allocated_reason(std::string* reason) {
  if (reason != nullptr) {
    _impl_._has_bits_[0] |= 0x00000001u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000001u;
  }
  _impl_.reason_.SetAllocated(reason, GetArenaForAllocation());
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.reason_.IsDefault()) {
    _impl_.reason_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  // @@protoc_insertion_point(field_set_allocated:apollo.common.EngageAdvice.reason)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace common
}  // namespace apollo

PROTOBUF_NAMESPACE_OPEN

template <> struct is_proto_enum< ::apollo::common::EngageAdvice_Advice> : ::std::true_type {};
template <>
inline const EnumDescriptor* GetEnumDescriptor< ::apollo::common::EngageAdvice_Advice>() {
  return ::apollo::common::EngageAdvice_Advice_descriptor();
}

PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_modules_2fcommon_2fproto_2fdrive_5fstate_2eproto
