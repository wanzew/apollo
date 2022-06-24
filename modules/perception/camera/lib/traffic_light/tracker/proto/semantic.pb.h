// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/perception/camera/lib/traffic_light/tracker/proto/semantic.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2ftracker_2fproto_2fsemantic_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2ftracker_2fproto_2fsemantic_2eproto

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
#define PROTOBUF_INTERNAL_EXPORT_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2ftracker_2fproto_2fsemantic_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2ftracker_2fproto_2fsemantic_2eproto {
  static const uint32_t offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2ftracker_2fproto_2fsemantic_2eproto;
namespace apollo {
namespace perception {
namespace camera {
namespace traffic_light {
namespace tracker {
class SemanticReviseParam;
struct SemanticReviseParamDefaultTypeInternal;
extern SemanticReviseParamDefaultTypeInternal _SemanticReviseParam_default_instance_;
}  // namespace tracker
}  // namespace traffic_light
}  // namespace camera
}  // namespace perception
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> ::apollo::perception::camera::traffic_light::tracker::SemanticReviseParam* Arena::CreateMaybeMessage<::apollo::perception::camera::traffic_light::tracker::SemanticReviseParam>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace apollo {
namespace perception {
namespace camera {
namespace traffic_light {
namespace tracker {

// ===================================================================

class SemanticReviseParam final :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:apollo.perception.camera.traffic_light.tracker.SemanticReviseParam) */ {
 public:
  inline SemanticReviseParam() : SemanticReviseParam(nullptr) {}
  ~SemanticReviseParam() override;
  explicit PROTOBUF_CONSTEXPR SemanticReviseParam(::PROTOBUF_NAMESPACE_ID::internal::ConstantInitialized);

  SemanticReviseParam(const SemanticReviseParam& from);
  SemanticReviseParam(SemanticReviseParam&& from) noexcept
    : SemanticReviseParam() {
    *this = ::std::move(from);
  }

  inline SemanticReviseParam& operator=(const SemanticReviseParam& from) {
    CopyFrom(from);
    return *this;
  }
  inline SemanticReviseParam& operator=(SemanticReviseParam&& from) noexcept {
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
  static const SemanticReviseParam& default_instance() {
    return *internal_default_instance();
  }
  static inline const SemanticReviseParam* internal_default_instance() {
    return reinterpret_cast<const SemanticReviseParam*>(
               &_SemanticReviseParam_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(SemanticReviseParam& a, SemanticReviseParam& b) {
    a.Swap(&b);
  }
  inline void Swap(SemanticReviseParam* other) {
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
  void UnsafeArenaSwap(SemanticReviseParam* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetOwningArena() == other->GetOwningArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  SemanticReviseParam* New(::PROTOBUF_NAMESPACE_ID::Arena* arena = nullptr) const final {
    return CreateMaybeMessage<SemanticReviseParam>(arena);
  }
  using ::PROTOBUF_NAMESPACE_ID::Message::CopyFrom;
  void CopyFrom(const SemanticReviseParam& from);
  using ::PROTOBUF_NAMESPACE_ID::Message::MergeFrom;
  void MergeFrom(const SemanticReviseParam& from);
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
  void InternalSwap(SemanticReviseParam* other);

  private:
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "apollo.perception.camera.traffic_light.tracker.SemanticReviseParam";
  }
  protected:
  explicit SemanticReviseParam(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                       bool is_message_owned = false);
  public:

  static const ClassData _class_data_;
  const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*GetClassData() const final;

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kHystereticThresholdCountFieldNumber = 3,
    kReviseTimeSecondFieldNumber = 1,
    kBlinkThresholdSecondFieldNumber = 2,
  };
  // optional int32 hysteretic_threshold_count = 3 [default = 2];
  bool has_hysteretic_threshold_count() const;
  private:
  bool _internal_has_hysteretic_threshold_count() const;
  public:
  void clear_hysteretic_threshold_count();
  int32_t hysteretic_threshold_count() const;
  void set_hysteretic_threshold_count(int32_t value);
  private:
  int32_t _internal_hysteretic_threshold_count() const;
  void _internal_set_hysteretic_threshold_count(int32_t value);
  public:

  // optional float revise_time_second = 1 [default = 1.5];
  bool has_revise_time_second() const;
  private:
  bool _internal_has_revise_time_second() const;
  public:
  void clear_revise_time_second();
  float revise_time_second() const;
  void set_revise_time_second(float value);
  private:
  float _internal_revise_time_second() const;
  void _internal_set_revise_time_second(float value);
  public:

  // optional float blink_threshold_second = 2 [default = 0.4];
  bool has_blink_threshold_second() const;
  private:
  bool _internal_has_blink_threshold_second() const;
  public:
  void clear_blink_threshold_second();
  float blink_threshold_second() const;
  void set_blink_threshold_second(float value);
  private:
  float _internal_blink_threshold_second() const;
  void _internal_set_blink_threshold_second(float value);
  public:

  // @@protoc_insertion_point(class_scope:apollo.perception.camera.traffic_light.tracker.SemanticReviseParam)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  struct Impl_ {
    ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
    mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
    int32_t hysteretic_threshold_count_;
    float revise_time_second_;
    float blink_threshold_second_;
  };
  union { Impl_ _impl_; };
  friend struct ::TableStruct_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2ftracker_2fproto_2fsemantic_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// SemanticReviseParam

// optional float revise_time_second = 1 [default = 1.5];
inline bool SemanticReviseParam::_internal_has_revise_time_second() const {
  bool value = (_impl_._has_bits_[0] & 0x00000002u) != 0;
  return value;
}
inline bool SemanticReviseParam::has_revise_time_second() const {
  return _internal_has_revise_time_second();
}
inline void SemanticReviseParam::clear_revise_time_second() {
  _impl_.revise_time_second_ = 1.5f;
  _impl_._has_bits_[0] &= ~0x00000002u;
}
inline float SemanticReviseParam::_internal_revise_time_second() const {
  return _impl_.revise_time_second_;
}
inline float SemanticReviseParam::revise_time_second() const {
  // @@protoc_insertion_point(field_get:apollo.perception.camera.traffic_light.tracker.SemanticReviseParam.revise_time_second)
  return _internal_revise_time_second();
}
inline void SemanticReviseParam::_internal_set_revise_time_second(float value) {
  _impl_._has_bits_[0] |= 0x00000002u;
  _impl_.revise_time_second_ = value;
}
inline void SemanticReviseParam::set_revise_time_second(float value) {
  _internal_set_revise_time_second(value);
  // @@protoc_insertion_point(field_set:apollo.perception.camera.traffic_light.tracker.SemanticReviseParam.revise_time_second)
}

// optional float blink_threshold_second = 2 [default = 0.4];
inline bool SemanticReviseParam::_internal_has_blink_threshold_second() const {
  bool value = (_impl_._has_bits_[0] & 0x00000004u) != 0;
  return value;
}
inline bool SemanticReviseParam::has_blink_threshold_second() const {
  return _internal_has_blink_threshold_second();
}
inline void SemanticReviseParam::clear_blink_threshold_second() {
  _impl_.blink_threshold_second_ = 0.4f;
  _impl_._has_bits_[0] &= ~0x00000004u;
}
inline float SemanticReviseParam::_internal_blink_threshold_second() const {
  return _impl_.blink_threshold_second_;
}
inline float SemanticReviseParam::blink_threshold_second() const {
  // @@protoc_insertion_point(field_get:apollo.perception.camera.traffic_light.tracker.SemanticReviseParam.blink_threshold_second)
  return _internal_blink_threshold_second();
}
inline void SemanticReviseParam::_internal_set_blink_threshold_second(float value) {
  _impl_._has_bits_[0] |= 0x00000004u;
  _impl_.blink_threshold_second_ = value;
}
inline void SemanticReviseParam::set_blink_threshold_second(float value) {
  _internal_set_blink_threshold_second(value);
  // @@protoc_insertion_point(field_set:apollo.perception.camera.traffic_light.tracker.SemanticReviseParam.blink_threshold_second)
}

// optional int32 hysteretic_threshold_count = 3 [default = 2];
inline bool SemanticReviseParam::_internal_has_hysteretic_threshold_count() const {
  bool value = (_impl_._has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool SemanticReviseParam::has_hysteretic_threshold_count() const {
  return _internal_has_hysteretic_threshold_count();
}
inline void SemanticReviseParam::clear_hysteretic_threshold_count() {
  _impl_.hysteretic_threshold_count_ = 2;
  _impl_._has_bits_[0] &= ~0x00000001u;
}
inline int32_t SemanticReviseParam::_internal_hysteretic_threshold_count() const {
  return _impl_.hysteretic_threshold_count_;
}
inline int32_t SemanticReviseParam::hysteretic_threshold_count() const {
  // @@protoc_insertion_point(field_get:apollo.perception.camera.traffic_light.tracker.SemanticReviseParam.hysteretic_threshold_count)
  return _internal_hysteretic_threshold_count();
}
inline void SemanticReviseParam::_internal_set_hysteretic_threshold_count(int32_t value) {
  _impl_._has_bits_[0] |= 0x00000001u;
  _impl_.hysteretic_threshold_count_ = value;
}
inline void SemanticReviseParam::set_hysteretic_threshold_count(int32_t value) {
  _internal_set_hysteretic_threshold_count(value);
  // @@protoc_insertion_point(field_set:apollo.perception.camera.traffic_light.tracker.SemanticReviseParam.hysteretic_threshold_count)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace tracker
}  // namespace traffic_light
}  // namespace camera
}  // namespace perception
}  // namespace apollo

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2ftracker_2fproto_2fsemantic_2eproto
