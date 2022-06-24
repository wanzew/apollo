// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/perception/onboard/proto/radar_component_config.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fonboard_2fproto_2fradar_5fcomponent_5fconfig_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fonboard_2fproto_2fradar_5fcomponent_5fconfig_2eproto

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
#define PROTOBUF_INTERNAL_EXPORT_modules_2fperception_2fonboard_2fproto_2fradar_5fcomponent_5fconfig_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_modules_2fperception_2fonboard_2fproto_2fradar_5fcomponent_5fconfig_2eproto {
  static const uint32_t offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_modules_2fperception_2fonboard_2fproto_2fradar_5fcomponent_5fconfig_2eproto;
namespace apollo {
namespace perception {
namespace onboard {
class RadarComponentConfig;
struct RadarComponentConfigDefaultTypeInternal;
extern RadarComponentConfigDefaultTypeInternal _RadarComponentConfig_default_instance_;
}  // namespace onboard
}  // namespace perception
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> ::apollo::perception::onboard::RadarComponentConfig* Arena::CreateMaybeMessage<::apollo::perception::onboard::RadarComponentConfig>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace apollo {
namespace perception {
namespace onboard {

// ===================================================================

class RadarComponentConfig final :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:apollo.perception.onboard.RadarComponentConfig) */ {
 public:
  inline RadarComponentConfig() : RadarComponentConfig(nullptr) {}
  ~RadarComponentConfig() override;
  explicit PROTOBUF_CONSTEXPR RadarComponentConfig(::PROTOBUF_NAMESPACE_ID::internal::ConstantInitialized);

  RadarComponentConfig(const RadarComponentConfig& from);
  RadarComponentConfig(RadarComponentConfig&& from) noexcept
    : RadarComponentConfig() {
    *this = ::std::move(from);
  }

  inline RadarComponentConfig& operator=(const RadarComponentConfig& from) {
    CopyFrom(from);
    return *this;
  }
  inline RadarComponentConfig& operator=(RadarComponentConfig&& from) noexcept {
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
  static const RadarComponentConfig& default_instance() {
    return *internal_default_instance();
  }
  static inline const RadarComponentConfig* internal_default_instance() {
    return reinterpret_cast<const RadarComponentConfig*>(
               &_RadarComponentConfig_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(RadarComponentConfig& a, RadarComponentConfig& b) {
    a.Swap(&b);
  }
  inline void Swap(RadarComponentConfig* other) {
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
  void UnsafeArenaSwap(RadarComponentConfig* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetOwningArena() == other->GetOwningArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  RadarComponentConfig* New(::PROTOBUF_NAMESPACE_ID::Arena* arena = nullptr) const final {
    return CreateMaybeMessage<RadarComponentConfig>(arena);
  }
  using ::PROTOBUF_NAMESPACE_ID::Message::CopyFrom;
  void CopyFrom(const RadarComponentConfig& from);
  using ::PROTOBUF_NAMESPACE_ID::Message::MergeFrom;
  void MergeFrom(const RadarComponentConfig& from);
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
  void InternalSwap(RadarComponentConfig* other);

  private:
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "apollo.perception.onboard.RadarComponentConfig";
  }
  protected:
  explicit RadarComponentConfig(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                       bool is_message_owned = false);
  public:

  static const ClassData _class_data_;
  const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*GetClassData() const final;

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kRadarNameFieldNumber = 1,
    kTfChildFrameIdFieldNumber = 2,
    kRadarPreprocessorMethodFieldNumber = 4,
    kRadarPerceptionMethodFieldNumber = 5,
    kRadarPipelineNameFieldNumber = 6,
    kOdometryChannelNameFieldNumber = 7,
    kOutputChannelNameFieldNumber = 8,
    kRadarForwardDistanceFieldNumber = 3,
  };
  // optional string radar_name = 1;
  bool has_radar_name() const;
  private:
  bool _internal_has_radar_name() const;
  public:
  void clear_radar_name();
  const std::string& radar_name() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_radar_name(ArgT0&& arg0, ArgT... args);
  std::string* mutable_radar_name();
  PROTOBUF_NODISCARD std::string* release_radar_name();
  void set_allocated_radar_name(std::string* radar_name);
  private:
  const std::string& _internal_radar_name() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_radar_name(const std::string& value);
  std::string* _internal_mutable_radar_name();
  public:

  // optional string tf_child_frame_id = 2;
  bool has_tf_child_frame_id() const;
  private:
  bool _internal_has_tf_child_frame_id() const;
  public:
  void clear_tf_child_frame_id();
  const std::string& tf_child_frame_id() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_tf_child_frame_id(ArgT0&& arg0, ArgT... args);
  std::string* mutable_tf_child_frame_id();
  PROTOBUF_NODISCARD std::string* release_tf_child_frame_id();
  void set_allocated_tf_child_frame_id(std::string* tf_child_frame_id);
  private:
  const std::string& _internal_tf_child_frame_id() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_tf_child_frame_id(const std::string& value);
  std::string* _internal_mutable_tf_child_frame_id();
  public:

  // optional string radar_preprocessor_method = 4;
  bool has_radar_preprocessor_method() const;
  private:
  bool _internal_has_radar_preprocessor_method() const;
  public:
  void clear_radar_preprocessor_method();
  const std::string& radar_preprocessor_method() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_radar_preprocessor_method(ArgT0&& arg0, ArgT... args);
  std::string* mutable_radar_preprocessor_method();
  PROTOBUF_NODISCARD std::string* release_radar_preprocessor_method();
  void set_allocated_radar_preprocessor_method(std::string* radar_preprocessor_method);
  private:
  const std::string& _internal_radar_preprocessor_method() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_radar_preprocessor_method(const std::string& value);
  std::string* _internal_mutable_radar_preprocessor_method();
  public:

  // optional string radar_perception_method = 5;
  bool has_radar_perception_method() const;
  private:
  bool _internal_has_radar_perception_method() const;
  public:
  void clear_radar_perception_method();
  const std::string& radar_perception_method() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_radar_perception_method(ArgT0&& arg0, ArgT... args);
  std::string* mutable_radar_perception_method();
  PROTOBUF_NODISCARD std::string* release_radar_perception_method();
  void set_allocated_radar_perception_method(std::string* radar_perception_method);
  private:
  const std::string& _internal_radar_perception_method() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_radar_perception_method(const std::string& value);
  std::string* _internal_mutable_radar_perception_method();
  public:

  // optional string radar_pipeline_name = 6;
  bool has_radar_pipeline_name() const;
  private:
  bool _internal_has_radar_pipeline_name() const;
  public:
  void clear_radar_pipeline_name();
  const std::string& radar_pipeline_name() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_radar_pipeline_name(ArgT0&& arg0, ArgT... args);
  std::string* mutable_radar_pipeline_name();
  PROTOBUF_NODISCARD std::string* release_radar_pipeline_name();
  void set_allocated_radar_pipeline_name(std::string* radar_pipeline_name);
  private:
  const std::string& _internal_radar_pipeline_name() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_radar_pipeline_name(const std::string& value);
  std::string* _internal_mutable_radar_pipeline_name();
  public:

  // optional string odometry_channel_name = 7;
  bool has_odometry_channel_name() const;
  private:
  bool _internal_has_odometry_channel_name() const;
  public:
  void clear_odometry_channel_name();
  const std::string& odometry_channel_name() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_odometry_channel_name(ArgT0&& arg0, ArgT... args);
  std::string* mutable_odometry_channel_name();
  PROTOBUF_NODISCARD std::string* release_odometry_channel_name();
  void set_allocated_odometry_channel_name(std::string* odometry_channel_name);
  private:
  const std::string& _internal_odometry_channel_name() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_odometry_channel_name(const std::string& value);
  std::string* _internal_mutable_odometry_channel_name();
  public:

  // optional string output_channel_name = 8;
  bool has_output_channel_name() const;
  private:
  bool _internal_has_output_channel_name() const;
  public:
  void clear_output_channel_name();
  const std::string& output_channel_name() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_output_channel_name(ArgT0&& arg0, ArgT... args);
  std::string* mutable_output_channel_name();
  PROTOBUF_NODISCARD std::string* release_output_channel_name();
  void set_allocated_output_channel_name(std::string* output_channel_name);
  private:
  const std::string& _internal_output_channel_name() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_output_channel_name(const std::string& value);
  std::string* _internal_mutable_output_channel_name();
  public:

  // optional double radar_forward_distance = 3;
  bool has_radar_forward_distance() const;
  private:
  bool _internal_has_radar_forward_distance() const;
  public:
  void clear_radar_forward_distance();
  double radar_forward_distance() const;
  void set_radar_forward_distance(double value);
  private:
  double _internal_radar_forward_distance() const;
  void _internal_set_radar_forward_distance(double value);
  public:

  // @@protoc_insertion_point(class_scope:apollo.perception.onboard.RadarComponentConfig)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  struct Impl_ {
    ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
    mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr radar_name_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr tf_child_frame_id_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr radar_preprocessor_method_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr radar_perception_method_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr radar_pipeline_name_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr odometry_channel_name_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr output_channel_name_;
    double radar_forward_distance_;
  };
  union { Impl_ _impl_; };
  friend struct ::TableStruct_modules_2fperception_2fonboard_2fproto_2fradar_5fcomponent_5fconfig_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// RadarComponentConfig

// optional string radar_name = 1;
inline bool RadarComponentConfig::_internal_has_radar_name() const {
  bool value = (_impl_._has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool RadarComponentConfig::has_radar_name() const {
  return _internal_has_radar_name();
}
inline void RadarComponentConfig::clear_radar_name() {
  _impl_.radar_name_.ClearToEmpty();
  _impl_._has_bits_[0] &= ~0x00000001u;
}
inline const std::string& RadarComponentConfig::radar_name() const {
  // @@protoc_insertion_point(field_get:apollo.perception.onboard.RadarComponentConfig.radar_name)
  return _internal_radar_name();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void RadarComponentConfig::set_radar_name(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000001u;
 _impl_.radar_name_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.onboard.RadarComponentConfig.radar_name)
}
inline std::string* RadarComponentConfig::mutable_radar_name() {
  std::string* _s = _internal_mutable_radar_name();
  // @@protoc_insertion_point(field_mutable:apollo.perception.onboard.RadarComponentConfig.radar_name)
  return _s;
}
inline const std::string& RadarComponentConfig::_internal_radar_name() const {
  return _impl_.radar_name_.Get();
}
inline void RadarComponentConfig::_internal_set_radar_name(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000001u;
  _impl_.radar_name_.Set(value, GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::_internal_mutable_radar_name() {
  _impl_._has_bits_[0] |= 0x00000001u;
  return _impl_.radar_name_.Mutable(GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::release_radar_name() {
  // @@protoc_insertion_point(field_release:apollo.perception.onboard.RadarComponentConfig.radar_name)
  if (!_internal_has_radar_name()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000001u;
  auto* p = _impl_.radar_name_.Release();
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.radar_name_.IsDefault()) {
    _impl_.radar_name_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  return p;
}
inline void RadarComponentConfig::set_allocated_radar_name(std::string* radar_name) {
  if (radar_name != nullptr) {
    _impl_._has_bits_[0] |= 0x00000001u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000001u;
  }
  _impl_.radar_name_.SetAllocated(radar_name, GetArenaForAllocation());
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.radar_name_.IsDefault()) {
    _impl_.radar_name_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.onboard.RadarComponentConfig.radar_name)
}

// optional string tf_child_frame_id = 2;
inline bool RadarComponentConfig::_internal_has_tf_child_frame_id() const {
  bool value = (_impl_._has_bits_[0] & 0x00000002u) != 0;
  return value;
}
inline bool RadarComponentConfig::has_tf_child_frame_id() const {
  return _internal_has_tf_child_frame_id();
}
inline void RadarComponentConfig::clear_tf_child_frame_id() {
  _impl_.tf_child_frame_id_.ClearToEmpty();
  _impl_._has_bits_[0] &= ~0x00000002u;
}
inline const std::string& RadarComponentConfig::tf_child_frame_id() const {
  // @@protoc_insertion_point(field_get:apollo.perception.onboard.RadarComponentConfig.tf_child_frame_id)
  return _internal_tf_child_frame_id();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void RadarComponentConfig::set_tf_child_frame_id(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000002u;
 _impl_.tf_child_frame_id_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.onboard.RadarComponentConfig.tf_child_frame_id)
}
inline std::string* RadarComponentConfig::mutable_tf_child_frame_id() {
  std::string* _s = _internal_mutable_tf_child_frame_id();
  // @@protoc_insertion_point(field_mutable:apollo.perception.onboard.RadarComponentConfig.tf_child_frame_id)
  return _s;
}
inline const std::string& RadarComponentConfig::_internal_tf_child_frame_id() const {
  return _impl_.tf_child_frame_id_.Get();
}
inline void RadarComponentConfig::_internal_set_tf_child_frame_id(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000002u;
  _impl_.tf_child_frame_id_.Set(value, GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::_internal_mutable_tf_child_frame_id() {
  _impl_._has_bits_[0] |= 0x00000002u;
  return _impl_.tf_child_frame_id_.Mutable(GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::release_tf_child_frame_id() {
  // @@protoc_insertion_point(field_release:apollo.perception.onboard.RadarComponentConfig.tf_child_frame_id)
  if (!_internal_has_tf_child_frame_id()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000002u;
  auto* p = _impl_.tf_child_frame_id_.Release();
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.tf_child_frame_id_.IsDefault()) {
    _impl_.tf_child_frame_id_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  return p;
}
inline void RadarComponentConfig::set_allocated_tf_child_frame_id(std::string* tf_child_frame_id) {
  if (tf_child_frame_id != nullptr) {
    _impl_._has_bits_[0] |= 0x00000002u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000002u;
  }
  _impl_.tf_child_frame_id_.SetAllocated(tf_child_frame_id, GetArenaForAllocation());
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.tf_child_frame_id_.IsDefault()) {
    _impl_.tf_child_frame_id_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.onboard.RadarComponentConfig.tf_child_frame_id)
}

// optional double radar_forward_distance = 3;
inline bool RadarComponentConfig::_internal_has_radar_forward_distance() const {
  bool value = (_impl_._has_bits_[0] & 0x00000080u) != 0;
  return value;
}
inline bool RadarComponentConfig::has_radar_forward_distance() const {
  return _internal_has_radar_forward_distance();
}
inline void RadarComponentConfig::clear_radar_forward_distance() {
  _impl_.radar_forward_distance_ = 0;
  _impl_._has_bits_[0] &= ~0x00000080u;
}
inline double RadarComponentConfig::_internal_radar_forward_distance() const {
  return _impl_.radar_forward_distance_;
}
inline double RadarComponentConfig::radar_forward_distance() const {
  // @@protoc_insertion_point(field_get:apollo.perception.onboard.RadarComponentConfig.radar_forward_distance)
  return _internal_radar_forward_distance();
}
inline void RadarComponentConfig::_internal_set_radar_forward_distance(double value) {
  _impl_._has_bits_[0] |= 0x00000080u;
  _impl_.radar_forward_distance_ = value;
}
inline void RadarComponentConfig::set_radar_forward_distance(double value) {
  _internal_set_radar_forward_distance(value);
  // @@protoc_insertion_point(field_set:apollo.perception.onboard.RadarComponentConfig.radar_forward_distance)
}

// optional string radar_preprocessor_method = 4;
inline bool RadarComponentConfig::_internal_has_radar_preprocessor_method() const {
  bool value = (_impl_._has_bits_[0] & 0x00000004u) != 0;
  return value;
}
inline bool RadarComponentConfig::has_radar_preprocessor_method() const {
  return _internal_has_radar_preprocessor_method();
}
inline void RadarComponentConfig::clear_radar_preprocessor_method() {
  _impl_.radar_preprocessor_method_.ClearToEmpty();
  _impl_._has_bits_[0] &= ~0x00000004u;
}
inline const std::string& RadarComponentConfig::radar_preprocessor_method() const {
  // @@protoc_insertion_point(field_get:apollo.perception.onboard.RadarComponentConfig.radar_preprocessor_method)
  return _internal_radar_preprocessor_method();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void RadarComponentConfig::set_radar_preprocessor_method(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000004u;
 _impl_.radar_preprocessor_method_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.onboard.RadarComponentConfig.radar_preprocessor_method)
}
inline std::string* RadarComponentConfig::mutable_radar_preprocessor_method() {
  std::string* _s = _internal_mutable_radar_preprocessor_method();
  // @@protoc_insertion_point(field_mutable:apollo.perception.onboard.RadarComponentConfig.radar_preprocessor_method)
  return _s;
}
inline const std::string& RadarComponentConfig::_internal_radar_preprocessor_method() const {
  return _impl_.radar_preprocessor_method_.Get();
}
inline void RadarComponentConfig::_internal_set_radar_preprocessor_method(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000004u;
  _impl_.radar_preprocessor_method_.Set(value, GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::_internal_mutable_radar_preprocessor_method() {
  _impl_._has_bits_[0] |= 0x00000004u;
  return _impl_.radar_preprocessor_method_.Mutable(GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::release_radar_preprocessor_method() {
  // @@protoc_insertion_point(field_release:apollo.perception.onboard.RadarComponentConfig.radar_preprocessor_method)
  if (!_internal_has_radar_preprocessor_method()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000004u;
  auto* p = _impl_.radar_preprocessor_method_.Release();
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.radar_preprocessor_method_.IsDefault()) {
    _impl_.radar_preprocessor_method_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  return p;
}
inline void RadarComponentConfig::set_allocated_radar_preprocessor_method(std::string* radar_preprocessor_method) {
  if (radar_preprocessor_method != nullptr) {
    _impl_._has_bits_[0] |= 0x00000004u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000004u;
  }
  _impl_.radar_preprocessor_method_.SetAllocated(radar_preprocessor_method, GetArenaForAllocation());
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.radar_preprocessor_method_.IsDefault()) {
    _impl_.radar_preprocessor_method_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.onboard.RadarComponentConfig.radar_preprocessor_method)
}

// optional string radar_perception_method = 5;
inline bool RadarComponentConfig::_internal_has_radar_perception_method() const {
  bool value = (_impl_._has_bits_[0] & 0x00000008u) != 0;
  return value;
}
inline bool RadarComponentConfig::has_radar_perception_method() const {
  return _internal_has_radar_perception_method();
}
inline void RadarComponentConfig::clear_radar_perception_method() {
  _impl_.radar_perception_method_.ClearToEmpty();
  _impl_._has_bits_[0] &= ~0x00000008u;
}
inline const std::string& RadarComponentConfig::radar_perception_method() const {
  // @@protoc_insertion_point(field_get:apollo.perception.onboard.RadarComponentConfig.radar_perception_method)
  return _internal_radar_perception_method();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void RadarComponentConfig::set_radar_perception_method(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000008u;
 _impl_.radar_perception_method_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.onboard.RadarComponentConfig.radar_perception_method)
}
inline std::string* RadarComponentConfig::mutable_radar_perception_method() {
  std::string* _s = _internal_mutable_radar_perception_method();
  // @@protoc_insertion_point(field_mutable:apollo.perception.onboard.RadarComponentConfig.radar_perception_method)
  return _s;
}
inline const std::string& RadarComponentConfig::_internal_radar_perception_method() const {
  return _impl_.radar_perception_method_.Get();
}
inline void RadarComponentConfig::_internal_set_radar_perception_method(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000008u;
  _impl_.radar_perception_method_.Set(value, GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::_internal_mutable_radar_perception_method() {
  _impl_._has_bits_[0] |= 0x00000008u;
  return _impl_.radar_perception_method_.Mutable(GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::release_radar_perception_method() {
  // @@protoc_insertion_point(field_release:apollo.perception.onboard.RadarComponentConfig.radar_perception_method)
  if (!_internal_has_radar_perception_method()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000008u;
  auto* p = _impl_.radar_perception_method_.Release();
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.radar_perception_method_.IsDefault()) {
    _impl_.radar_perception_method_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  return p;
}
inline void RadarComponentConfig::set_allocated_radar_perception_method(std::string* radar_perception_method) {
  if (radar_perception_method != nullptr) {
    _impl_._has_bits_[0] |= 0x00000008u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000008u;
  }
  _impl_.radar_perception_method_.SetAllocated(radar_perception_method, GetArenaForAllocation());
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.radar_perception_method_.IsDefault()) {
    _impl_.radar_perception_method_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.onboard.RadarComponentConfig.radar_perception_method)
}

// optional string radar_pipeline_name = 6;
inline bool RadarComponentConfig::_internal_has_radar_pipeline_name() const {
  bool value = (_impl_._has_bits_[0] & 0x00000010u) != 0;
  return value;
}
inline bool RadarComponentConfig::has_radar_pipeline_name() const {
  return _internal_has_radar_pipeline_name();
}
inline void RadarComponentConfig::clear_radar_pipeline_name() {
  _impl_.radar_pipeline_name_.ClearToEmpty();
  _impl_._has_bits_[0] &= ~0x00000010u;
}
inline const std::string& RadarComponentConfig::radar_pipeline_name() const {
  // @@protoc_insertion_point(field_get:apollo.perception.onboard.RadarComponentConfig.radar_pipeline_name)
  return _internal_radar_pipeline_name();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void RadarComponentConfig::set_radar_pipeline_name(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000010u;
 _impl_.radar_pipeline_name_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.onboard.RadarComponentConfig.radar_pipeline_name)
}
inline std::string* RadarComponentConfig::mutable_radar_pipeline_name() {
  std::string* _s = _internal_mutable_radar_pipeline_name();
  // @@protoc_insertion_point(field_mutable:apollo.perception.onboard.RadarComponentConfig.radar_pipeline_name)
  return _s;
}
inline const std::string& RadarComponentConfig::_internal_radar_pipeline_name() const {
  return _impl_.radar_pipeline_name_.Get();
}
inline void RadarComponentConfig::_internal_set_radar_pipeline_name(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000010u;
  _impl_.radar_pipeline_name_.Set(value, GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::_internal_mutable_radar_pipeline_name() {
  _impl_._has_bits_[0] |= 0x00000010u;
  return _impl_.radar_pipeline_name_.Mutable(GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::release_radar_pipeline_name() {
  // @@protoc_insertion_point(field_release:apollo.perception.onboard.RadarComponentConfig.radar_pipeline_name)
  if (!_internal_has_radar_pipeline_name()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000010u;
  auto* p = _impl_.radar_pipeline_name_.Release();
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.radar_pipeline_name_.IsDefault()) {
    _impl_.radar_pipeline_name_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  return p;
}
inline void RadarComponentConfig::set_allocated_radar_pipeline_name(std::string* radar_pipeline_name) {
  if (radar_pipeline_name != nullptr) {
    _impl_._has_bits_[0] |= 0x00000010u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000010u;
  }
  _impl_.radar_pipeline_name_.SetAllocated(radar_pipeline_name, GetArenaForAllocation());
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.radar_pipeline_name_.IsDefault()) {
    _impl_.radar_pipeline_name_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.onboard.RadarComponentConfig.radar_pipeline_name)
}

// optional string odometry_channel_name = 7;
inline bool RadarComponentConfig::_internal_has_odometry_channel_name() const {
  bool value = (_impl_._has_bits_[0] & 0x00000020u) != 0;
  return value;
}
inline bool RadarComponentConfig::has_odometry_channel_name() const {
  return _internal_has_odometry_channel_name();
}
inline void RadarComponentConfig::clear_odometry_channel_name() {
  _impl_.odometry_channel_name_.ClearToEmpty();
  _impl_._has_bits_[0] &= ~0x00000020u;
}
inline const std::string& RadarComponentConfig::odometry_channel_name() const {
  // @@protoc_insertion_point(field_get:apollo.perception.onboard.RadarComponentConfig.odometry_channel_name)
  return _internal_odometry_channel_name();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void RadarComponentConfig::set_odometry_channel_name(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000020u;
 _impl_.odometry_channel_name_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.onboard.RadarComponentConfig.odometry_channel_name)
}
inline std::string* RadarComponentConfig::mutable_odometry_channel_name() {
  std::string* _s = _internal_mutable_odometry_channel_name();
  // @@protoc_insertion_point(field_mutable:apollo.perception.onboard.RadarComponentConfig.odometry_channel_name)
  return _s;
}
inline const std::string& RadarComponentConfig::_internal_odometry_channel_name() const {
  return _impl_.odometry_channel_name_.Get();
}
inline void RadarComponentConfig::_internal_set_odometry_channel_name(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000020u;
  _impl_.odometry_channel_name_.Set(value, GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::_internal_mutable_odometry_channel_name() {
  _impl_._has_bits_[0] |= 0x00000020u;
  return _impl_.odometry_channel_name_.Mutable(GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::release_odometry_channel_name() {
  // @@protoc_insertion_point(field_release:apollo.perception.onboard.RadarComponentConfig.odometry_channel_name)
  if (!_internal_has_odometry_channel_name()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000020u;
  auto* p = _impl_.odometry_channel_name_.Release();
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.odometry_channel_name_.IsDefault()) {
    _impl_.odometry_channel_name_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  return p;
}
inline void RadarComponentConfig::set_allocated_odometry_channel_name(std::string* odometry_channel_name) {
  if (odometry_channel_name != nullptr) {
    _impl_._has_bits_[0] |= 0x00000020u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000020u;
  }
  _impl_.odometry_channel_name_.SetAllocated(odometry_channel_name, GetArenaForAllocation());
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.odometry_channel_name_.IsDefault()) {
    _impl_.odometry_channel_name_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.onboard.RadarComponentConfig.odometry_channel_name)
}

// optional string output_channel_name = 8;
inline bool RadarComponentConfig::_internal_has_output_channel_name() const {
  bool value = (_impl_._has_bits_[0] & 0x00000040u) != 0;
  return value;
}
inline bool RadarComponentConfig::has_output_channel_name() const {
  return _internal_has_output_channel_name();
}
inline void RadarComponentConfig::clear_output_channel_name() {
  _impl_.output_channel_name_.ClearToEmpty();
  _impl_._has_bits_[0] &= ~0x00000040u;
}
inline const std::string& RadarComponentConfig::output_channel_name() const {
  // @@protoc_insertion_point(field_get:apollo.perception.onboard.RadarComponentConfig.output_channel_name)
  return _internal_output_channel_name();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void RadarComponentConfig::set_output_channel_name(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000040u;
 _impl_.output_channel_name_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.onboard.RadarComponentConfig.output_channel_name)
}
inline std::string* RadarComponentConfig::mutable_output_channel_name() {
  std::string* _s = _internal_mutable_output_channel_name();
  // @@protoc_insertion_point(field_mutable:apollo.perception.onboard.RadarComponentConfig.output_channel_name)
  return _s;
}
inline const std::string& RadarComponentConfig::_internal_output_channel_name() const {
  return _impl_.output_channel_name_.Get();
}
inline void RadarComponentConfig::_internal_set_output_channel_name(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000040u;
  _impl_.output_channel_name_.Set(value, GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::_internal_mutable_output_channel_name() {
  _impl_._has_bits_[0] |= 0x00000040u;
  return _impl_.output_channel_name_.Mutable(GetArenaForAllocation());
}
inline std::string* RadarComponentConfig::release_output_channel_name() {
  // @@protoc_insertion_point(field_release:apollo.perception.onboard.RadarComponentConfig.output_channel_name)
  if (!_internal_has_output_channel_name()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000040u;
  auto* p = _impl_.output_channel_name_.Release();
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.output_channel_name_.IsDefault()) {
    _impl_.output_channel_name_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  return p;
}
inline void RadarComponentConfig::set_allocated_output_channel_name(std::string* output_channel_name) {
  if (output_channel_name != nullptr) {
    _impl_._has_bits_[0] |= 0x00000040u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000040u;
  }
  _impl_.output_channel_name_.SetAllocated(output_channel_name, GetArenaForAllocation());
#ifdef PROTOBUF_FORCE_COPY_DEFAULT_STRING
  if (_impl_.output_channel_name_.IsDefault()) {
    _impl_.output_channel_name_.Set("", GetArenaForAllocation());
  }
#endif // PROTOBUF_FORCE_COPY_DEFAULT_STRING
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.onboard.RadarComponentConfig.output_channel_name)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace onboard
}  // namespace perception
}  // namespace apollo

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fonboard_2fproto_2fradar_5fcomponent_5fconfig_2eproto
