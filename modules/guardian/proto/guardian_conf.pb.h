// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/guardian/proto/guardian_conf.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_modules_2fguardian_2fproto_2fguardian_5fconf_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_modules_2fguardian_2fproto_2fguardian_5fconf_2eproto

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
#include <google/protobuf/unknown_field_set.h>
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_modules_2fguardian_2fproto_2fguardian_5fconf_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_modules_2fguardian_2fproto_2fguardian_5fconf_2eproto {
  static const uint32_t offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_modules_2fguardian_2fproto_2fguardian_5fconf_2eproto;
namespace apollo {
namespace guardian {
class GuardianConf;
struct GuardianConfDefaultTypeInternal;
extern GuardianConfDefaultTypeInternal _GuardianConf_default_instance_;
}  // namespace guardian
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> ::apollo::guardian::GuardianConf* Arena::CreateMaybeMessage<::apollo::guardian::GuardianConf>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace apollo {
namespace guardian {

// ===================================================================

class GuardianConf final :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:apollo.guardian.GuardianConf) */ {
 public:
  inline GuardianConf() : GuardianConf(nullptr) {}
  ~GuardianConf() override;
  explicit PROTOBUF_CONSTEXPR GuardianConf(::PROTOBUF_NAMESPACE_ID::internal::ConstantInitialized);

  GuardianConf(const GuardianConf& from);
  GuardianConf(GuardianConf&& from) noexcept
    : GuardianConf() {
    *this = ::std::move(from);
  }

  inline GuardianConf& operator=(const GuardianConf& from) {
    CopyFrom(from);
    return *this;
  }
  inline GuardianConf& operator=(GuardianConf&& from) noexcept {
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
  static const GuardianConf& default_instance() {
    return *internal_default_instance();
  }
  static inline const GuardianConf* internal_default_instance() {
    return reinterpret_cast<const GuardianConf*>(
               &_GuardianConf_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(GuardianConf& a, GuardianConf& b) {
    a.Swap(&b);
  }
  inline void Swap(GuardianConf* other) {
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
  void UnsafeArenaSwap(GuardianConf* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetOwningArena() == other->GetOwningArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  GuardianConf* New(::PROTOBUF_NAMESPACE_ID::Arena* arena = nullptr) const final {
    return CreateMaybeMessage<GuardianConf>(arena);
  }
  using ::PROTOBUF_NAMESPACE_ID::Message::CopyFrom;
  void CopyFrom(const GuardianConf& from);
  using ::PROTOBUF_NAMESPACE_ID::Message::MergeFrom;
  void MergeFrom(const GuardianConf& from);
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
  void InternalSwap(GuardianConf* other);

  private:
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "apollo.guardian.GuardianConf";
  }
  protected:
  explicit GuardianConf(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                       bool is_message_owned = false);
  public:

  static const ClassData _class_data_;
  const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*GetClassData() const final;

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kGuardianEnableFieldNumber = 1,
    kGuardianCmdEmergencyStopPercentageFieldNumber = 2,
    kGuardianCmdSoftStopPercentageFieldNumber = 3,
  };
  // optional bool guardian_enable = 1 [default = false];
  bool has_guardian_enable() const;
  private:
  bool _internal_has_guardian_enable() const;
  public:
  void clear_guardian_enable();
  bool guardian_enable() const;
  void set_guardian_enable(bool value);
  private:
  bool _internal_guardian_enable() const;
  void _internal_set_guardian_enable(bool value);
  public:

  // optional double guardian_cmd_emergency_stop_percentage = 2 [default = 50];
  bool has_guardian_cmd_emergency_stop_percentage() const;
  private:
  bool _internal_has_guardian_cmd_emergency_stop_percentage() const;
  public:
  void clear_guardian_cmd_emergency_stop_percentage();
  double guardian_cmd_emergency_stop_percentage() const;
  void set_guardian_cmd_emergency_stop_percentage(double value);
  private:
  double _internal_guardian_cmd_emergency_stop_percentage() const;
  void _internal_set_guardian_cmd_emergency_stop_percentage(double value);
  public:

  // optional double guardian_cmd_soft_stop_percentage = 3 [default = 25];
  bool has_guardian_cmd_soft_stop_percentage() const;
  private:
  bool _internal_has_guardian_cmd_soft_stop_percentage() const;
  public:
  void clear_guardian_cmd_soft_stop_percentage();
  double guardian_cmd_soft_stop_percentage() const;
  void set_guardian_cmd_soft_stop_percentage(double value);
  private:
  double _internal_guardian_cmd_soft_stop_percentage() const;
  void _internal_set_guardian_cmd_soft_stop_percentage(double value);
  public:

  // @@protoc_insertion_point(class_scope:apollo.guardian.GuardianConf)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  struct Impl_ {
    ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
    mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
    bool guardian_enable_;
    double guardian_cmd_emergency_stop_percentage_;
    double guardian_cmd_soft_stop_percentage_;
  };
  union { Impl_ _impl_; };
  friend struct ::TableStruct_modules_2fguardian_2fproto_2fguardian_5fconf_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// GuardianConf

// optional bool guardian_enable = 1 [default = false];
inline bool GuardianConf::_internal_has_guardian_enable() const {
  bool value = (_impl_._has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool GuardianConf::has_guardian_enable() const {
  return _internal_has_guardian_enable();
}
inline void GuardianConf::clear_guardian_enable() {
  _impl_.guardian_enable_ = false;
  _impl_._has_bits_[0] &= ~0x00000001u;
}
inline bool GuardianConf::_internal_guardian_enable() const {
  return _impl_.guardian_enable_;
}
inline bool GuardianConf::guardian_enable() const {
  // @@protoc_insertion_point(field_get:apollo.guardian.GuardianConf.guardian_enable)
  return _internal_guardian_enable();
}
inline void GuardianConf::_internal_set_guardian_enable(bool value) {
  _impl_._has_bits_[0] |= 0x00000001u;
  _impl_.guardian_enable_ = value;
}
inline void GuardianConf::set_guardian_enable(bool value) {
  _internal_set_guardian_enable(value);
  // @@protoc_insertion_point(field_set:apollo.guardian.GuardianConf.guardian_enable)
}

// optional double guardian_cmd_emergency_stop_percentage = 2 [default = 50];
inline bool GuardianConf::_internal_has_guardian_cmd_emergency_stop_percentage() const {
  bool value = (_impl_._has_bits_[0] & 0x00000002u) != 0;
  return value;
}
inline bool GuardianConf::has_guardian_cmd_emergency_stop_percentage() const {
  return _internal_has_guardian_cmd_emergency_stop_percentage();
}
inline void GuardianConf::clear_guardian_cmd_emergency_stop_percentage() {
  _impl_.guardian_cmd_emergency_stop_percentage_ = 50;
  _impl_._has_bits_[0] &= ~0x00000002u;
}
inline double GuardianConf::_internal_guardian_cmd_emergency_stop_percentage() const {
  return _impl_.guardian_cmd_emergency_stop_percentage_;
}
inline double GuardianConf::guardian_cmd_emergency_stop_percentage() const {
  // @@protoc_insertion_point(field_get:apollo.guardian.GuardianConf.guardian_cmd_emergency_stop_percentage)
  return _internal_guardian_cmd_emergency_stop_percentage();
}
inline void GuardianConf::_internal_set_guardian_cmd_emergency_stop_percentage(double value) {
  _impl_._has_bits_[0] |= 0x00000002u;
  _impl_.guardian_cmd_emergency_stop_percentage_ = value;
}
inline void GuardianConf::set_guardian_cmd_emergency_stop_percentage(double value) {
  _internal_set_guardian_cmd_emergency_stop_percentage(value);
  // @@protoc_insertion_point(field_set:apollo.guardian.GuardianConf.guardian_cmd_emergency_stop_percentage)
}

// optional double guardian_cmd_soft_stop_percentage = 3 [default = 25];
inline bool GuardianConf::_internal_has_guardian_cmd_soft_stop_percentage() const {
  bool value = (_impl_._has_bits_[0] & 0x00000004u) != 0;
  return value;
}
inline bool GuardianConf::has_guardian_cmd_soft_stop_percentage() const {
  return _internal_has_guardian_cmd_soft_stop_percentage();
}
inline void GuardianConf::clear_guardian_cmd_soft_stop_percentage() {
  _impl_.guardian_cmd_soft_stop_percentage_ = 25;
  _impl_._has_bits_[0] &= ~0x00000004u;
}
inline double GuardianConf::_internal_guardian_cmd_soft_stop_percentage() const {
  return _impl_.guardian_cmd_soft_stop_percentage_;
}
inline double GuardianConf::guardian_cmd_soft_stop_percentage() const {
  // @@protoc_insertion_point(field_get:apollo.guardian.GuardianConf.guardian_cmd_soft_stop_percentage)
  return _internal_guardian_cmd_soft_stop_percentage();
}
inline void GuardianConf::_internal_set_guardian_cmd_soft_stop_percentage(double value) {
  _impl_._has_bits_[0] |= 0x00000004u;
  _impl_.guardian_cmd_soft_stop_percentage_ = value;
}
inline void GuardianConf::set_guardian_cmd_soft_stop_percentage(double value) {
  _internal_set_guardian_cmd_soft_stop_percentage(value);
  // @@protoc_insertion_point(field_set:apollo.guardian.GuardianConf.guardian_cmd_soft_stop_percentage)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace guardian
}  // namespace apollo

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_modules_2fguardian_2fproto_2fguardian_5fconf_2eproto
