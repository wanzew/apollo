// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/perception/proto/async_fusion_config.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fproto_2fasync_5ffusion_5fconfig_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fproto_2fasync_5ffusion_5fconfig_2eproto

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
#define PROTOBUF_INTERNAL_EXPORT_modules_2fperception_2fproto_2fasync_5ffusion_5fconfig_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_modules_2fperception_2fproto_2fasync_5ffusion_5fconfig_2eproto {
  static const uint32_t offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_modules_2fperception_2fproto_2fasync_5ffusion_5fconfig_2eproto;
namespace apollo {
namespace perception {
namespace async_fusion_config {
class ModelConfigs;
struct ModelConfigsDefaultTypeInternal;
extern ModelConfigsDefaultTypeInternal _ModelConfigs_default_instance_;
}  // namespace async_fusion_config
}  // namespace perception
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> ::apollo::perception::async_fusion_config::ModelConfigs* Arena::CreateMaybeMessage<::apollo::perception::async_fusion_config::ModelConfigs>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace apollo {
namespace perception {
namespace async_fusion_config {

// ===================================================================

class ModelConfigs final :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:apollo.perception.async_fusion_config.ModelConfigs) */ {
 public:
  inline ModelConfigs() : ModelConfigs(nullptr) {}
  ~ModelConfigs() override;
  explicit PROTOBUF_CONSTEXPR ModelConfigs(::PROTOBUF_NAMESPACE_ID::internal::ConstantInitialized);

  ModelConfigs(const ModelConfigs& from);
  ModelConfigs(ModelConfigs&& from) noexcept
    : ModelConfigs() {
    *this = ::std::move(from);
  }

  inline ModelConfigs& operator=(const ModelConfigs& from) {
    CopyFrom(from);
    return *this;
  }
  inline ModelConfigs& operator=(ModelConfigs&& from) noexcept {
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
  static const ModelConfigs& default_instance() {
    return *internal_default_instance();
  }
  static inline const ModelConfigs* internal_default_instance() {
    return reinterpret_cast<const ModelConfigs*>(
               &_ModelConfigs_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(ModelConfigs& a, ModelConfigs& b) {
    a.Swap(&b);
  }
  inline void Swap(ModelConfigs* other) {
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
  void UnsafeArenaSwap(ModelConfigs* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetOwningArena() == other->GetOwningArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  ModelConfigs* New(::PROTOBUF_NAMESPACE_ID::Arena* arena = nullptr) const final {
    return CreateMaybeMessage<ModelConfigs>(arena);
  }
  using ::PROTOBUF_NAMESPACE_ID::Message::CopyFrom;
  void CopyFrom(const ModelConfigs& from);
  using ::PROTOBUF_NAMESPACE_ID::Message::MergeFrom;
  void MergeFrom(const ModelConfigs& from);
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
  void InternalSwap(ModelConfigs* other);

  private:
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "apollo.perception.async_fusion_config.ModelConfigs";
  }
  protected:
  explicit ModelConfigs(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                       bool is_message_owned = false);
  public:

  static const ClassData _class_data_;
  const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*GetClassData() const final;

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kNameFieldNumber = 1,
    kVersionFieldNumber = 2,
    kMatchMethodFieldNumber = 3,
    kPublishSensorFieldNumber = 11,
    kMaxMatchDistanceFieldNumber = 4,
    kMaxLidarInvisiblePeriodFieldNumber = 5,
    kMaxRadarInvisiblePeriodFieldNumber = 6,
    kMaxRadarConfidentAngleFieldNumber = 7,
    kMinRadarConfidentDistanceFieldNumber = 8,
    kPublishIfHasLidarFieldNumber = 9,
    kPublishIfHasRadarFieldNumber = 10,
    kUseRadarFieldNumber = 12,
    kUseLidarFieldNumber = 13,
  };
  // optional string name = 1 [default = "AsyncFusion"];
  bool has_name() const;
  private:
  bool _internal_has_name() const;
  public:
  void clear_name();
  const std::string& name() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_name(ArgT0&& arg0, ArgT... args);
  std::string* mutable_name();
  PROTOBUF_NODISCARD std::string* release_name();
  void set_allocated_name(std::string* name);
  private:
  const std::string& _internal_name() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_name(const std::string& value);
  std::string* _internal_mutable_name();
  public:

  // optional string version = 2 [default = "1.0.0"];
  bool has_version() const;
  private:
  bool _internal_has_version() const;
  public:
  void clear_version();
  const std::string& version() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_version(ArgT0&& arg0, ArgT... args);
  std::string* mutable_version();
  PROTOBUF_NODISCARD std::string* release_version();
  void set_allocated_version(std::string* version);
  private:
  const std::string& _internal_version() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_version(const std::string& value);
  std::string* _internal_mutable_version();
  public:

  // optional string match_method = 3 [default = "hm_matcher"];
  bool has_match_method() const;
  private:
  bool _internal_has_match_method() const;
  public:
  void clear_match_method();
  const std::string& match_method() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_match_method(ArgT0&& arg0, ArgT... args);
  std::string* mutable_match_method();
  PROTOBUF_NODISCARD std::string* release_match_method();
  void set_allocated_match_method(std::string* match_method);
  private:
  const std::string& _internal_match_method() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_match_method(const std::string& value);
  std::string* _internal_mutable_match_method();
  public:

  // optional string publish_sensor = 11 [default = "velodyne_64"];
  bool has_publish_sensor() const;
  private:
  bool _internal_has_publish_sensor() const;
  public:
  void clear_publish_sensor();
  const std::string& publish_sensor() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_publish_sensor(ArgT0&& arg0, ArgT... args);
  std::string* mutable_publish_sensor();
  PROTOBUF_NODISCARD std::string* release_publish_sensor();
  void set_allocated_publish_sensor(std::string* publish_sensor);
  private:
  const std::string& _internal_publish_sensor() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_publish_sensor(const std::string& value);
  std::string* _internal_mutable_publish_sensor();
  public:

  // optional float max_match_distance = 4 [default = 4];
  bool has_max_match_distance() const;
  private:
  bool _internal_has_max_match_distance() const;
  public:
  void clear_max_match_distance();
  float max_match_distance() const;
  void set_max_match_distance(float value);
  private:
  float _internal_max_match_distance() const;
  void _internal_set_max_match_distance(float value);
  public:

  // optional float max_lidar_invisible_period = 5 [default = 0.25];
  bool has_max_lidar_invisible_period() const;
  private:
  bool _internal_has_max_lidar_invisible_period() const;
  public:
  void clear_max_lidar_invisible_period();
  float max_lidar_invisible_period() const;
  void set_max_lidar_invisible_period(float value);
  private:
  float _internal_max_lidar_invisible_period() const;
  void _internal_set_max_lidar_invisible_period(float value);
  public:

  // optional float max_radar_invisible_period = 6 [default = 0.25];
  bool has_max_radar_invisible_period() const;
  private:
  bool _internal_has_max_radar_invisible_period() const;
  public:
  void clear_max_radar_invisible_period();
  float max_radar_invisible_period() const;
  void set_max_radar_invisible_period(float value);
  private:
  float _internal_max_radar_invisible_period() const;
  void _internal_set_max_radar_invisible_period(float value);
  public:

  // optional float max_radar_confident_angle = 7 [default = 30];
  bool has_max_radar_confident_angle() const;
  private:
  bool _internal_has_max_radar_confident_angle() const;
  public:
  void clear_max_radar_confident_angle();
  float max_radar_confident_angle() const;
  void set_max_radar_confident_angle(float value);
  private:
  float _internal_max_radar_confident_angle() const;
  void _internal_set_max_radar_confident_angle(float value);
  public:

  // optional float min_radar_confident_distance = 8 [default = 40];
  bool has_min_radar_confident_distance() const;
  private:
  bool _internal_has_min_radar_confident_distance() const;
  public:
  void clear_min_radar_confident_distance();
  float min_radar_confident_distance() const;
  void set_min_radar_confident_distance(float value);
  private:
  float _internal_min_radar_confident_distance() const;
  void _internal_set_min_radar_confident_distance(float value);
  public:

  // optional bool publish_if_has_lidar = 9 [default = true];
  bool has_publish_if_has_lidar() const;
  private:
  bool _internal_has_publish_if_has_lidar() const;
  public:
  void clear_publish_if_has_lidar();
  bool publish_if_has_lidar() const;
  void set_publish_if_has_lidar(bool value);
  private:
  bool _internal_publish_if_has_lidar() const;
  void _internal_set_publish_if_has_lidar(bool value);
  public:

  // optional bool publish_if_has_radar = 10 [default = true];
  bool has_publish_if_has_radar() const;
  private:
  bool _internal_has_publish_if_has_radar() const;
  public:
  void clear_publish_if_has_radar();
  bool publish_if_has_radar() const;
  void set_publish_if_has_radar(bool value);
  private:
  bool _internal_publish_if_has_radar() const;
  void _internal_set_publish_if_has_radar(bool value);
  public:

  // optional bool use_radar = 12 [default = true];
  bool has_use_radar() const;
  private:
  bool _internal_has_use_radar() const;
  public:
  void clear_use_radar();
  bool use_radar() const;
  void set_use_radar(bool value);
  private:
  bool _internal_use_radar() const;
  void _internal_set_use_radar(bool value);
  public:

  // optional bool use_lidar = 13 [default = true];
  bool has_use_lidar() const;
  private:
  bool _internal_has_use_lidar() const;
  public:
  void clear_use_lidar();
  bool use_lidar() const;
  void set_use_lidar(bool value);
  private:
  bool _internal_use_lidar() const;
  void _internal_set_use_lidar(bool value);
  public:

  // @@protoc_insertion_point(class_scope:apollo.perception.async_fusion_config.ModelConfigs)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  struct Impl_ {
    ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
    mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
    static const ::PROTOBUF_NAMESPACE_ID::internal::LazyString _i_give_permission_to_break_this_code_default_name_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr name_;
    static const ::PROTOBUF_NAMESPACE_ID::internal::LazyString _i_give_permission_to_break_this_code_default_version_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr version_;
    static const ::PROTOBUF_NAMESPACE_ID::internal::LazyString _i_give_permission_to_break_this_code_default_match_method_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr match_method_;
    static const ::PROTOBUF_NAMESPACE_ID::internal::LazyString _i_give_permission_to_break_this_code_default_publish_sensor_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr publish_sensor_;
    float max_match_distance_;
    float max_lidar_invisible_period_;
    float max_radar_invisible_period_;
    float max_radar_confident_angle_;
    float min_radar_confident_distance_;
    bool publish_if_has_lidar_;
    bool publish_if_has_radar_;
    bool use_radar_;
    bool use_lidar_;
  };
  union { Impl_ _impl_; };
  friend struct ::TableStruct_modules_2fperception_2fproto_2fasync_5ffusion_5fconfig_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// ModelConfigs

// optional string name = 1 [default = "AsyncFusion"];
inline bool ModelConfigs::_internal_has_name() const {
  bool value = (_impl_._has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool ModelConfigs::has_name() const {
  return _internal_has_name();
}
inline void ModelConfigs::clear_name() {
  _impl_.name_.ClearToDefault(::apollo::perception::async_fusion_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_name_, GetArenaForAllocation());
  _impl_._has_bits_[0] &= ~0x00000001u;
}
inline const std::string& ModelConfigs::name() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.name)
  if (_impl_.name_.IsDefault()) return Impl_::_i_give_permission_to_break_this_code_default_name_.get();
  return _internal_name();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void ModelConfigs::set_name(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000001u;
 _impl_.name_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.name)
}
inline std::string* ModelConfigs::mutable_name() {
  std::string* _s = _internal_mutable_name();
  // @@protoc_insertion_point(field_mutable:apollo.perception.async_fusion_config.ModelConfigs.name)
  return _s;
}
inline const std::string& ModelConfigs::_internal_name() const {
  return _impl_.name_.Get();
}
inline void ModelConfigs::_internal_set_name(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000001u;
  _impl_.name_.Set(value, GetArenaForAllocation());
}
inline std::string* ModelConfigs::_internal_mutable_name() {
  _impl_._has_bits_[0] |= 0x00000001u;
  return _impl_.name_.Mutable(::apollo::perception::async_fusion_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_name_, GetArenaForAllocation());
}
inline std::string* ModelConfigs::release_name() {
  // @@protoc_insertion_point(field_release:apollo.perception.async_fusion_config.ModelConfigs.name)
  if (!_internal_has_name()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000001u;
  auto* p = _impl_.name_.Release();
  return p;
}
inline void ModelConfigs::set_allocated_name(std::string* name) {
  if (name != nullptr) {
    _impl_._has_bits_[0] |= 0x00000001u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000001u;
  }
  _impl_.name_.SetAllocated(name, GetArenaForAllocation());
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.async_fusion_config.ModelConfigs.name)
}

// optional string version = 2 [default = "1.0.0"];
inline bool ModelConfigs::_internal_has_version() const {
  bool value = (_impl_._has_bits_[0] & 0x00000002u) != 0;
  return value;
}
inline bool ModelConfigs::has_version() const {
  return _internal_has_version();
}
inline void ModelConfigs::clear_version() {
  _impl_.version_.ClearToDefault(::apollo::perception::async_fusion_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_version_, GetArenaForAllocation());
  _impl_._has_bits_[0] &= ~0x00000002u;
}
inline const std::string& ModelConfigs::version() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.version)
  if (_impl_.version_.IsDefault()) return Impl_::_i_give_permission_to_break_this_code_default_version_.get();
  return _internal_version();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void ModelConfigs::set_version(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000002u;
 _impl_.version_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.version)
}
inline std::string* ModelConfigs::mutable_version() {
  std::string* _s = _internal_mutable_version();
  // @@protoc_insertion_point(field_mutable:apollo.perception.async_fusion_config.ModelConfigs.version)
  return _s;
}
inline const std::string& ModelConfigs::_internal_version() const {
  return _impl_.version_.Get();
}
inline void ModelConfigs::_internal_set_version(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000002u;
  _impl_.version_.Set(value, GetArenaForAllocation());
}
inline std::string* ModelConfigs::_internal_mutable_version() {
  _impl_._has_bits_[0] |= 0x00000002u;
  return _impl_.version_.Mutable(::apollo::perception::async_fusion_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_version_, GetArenaForAllocation());
}
inline std::string* ModelConfigs::release_version() {
  // @@protoc_insertion_point(field_release:apollo.perception.async_fusion_config.ModelConfigs.version)
  if (!_internal_has_version()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000002u;
  auto* p = _impl_.version_.Release();
  return p;
}
inline void ModelConfigs::set_allocated_version(std::string* version) {
  if (version != nullptr) {
    _impl_._has_bits_[0] |= 0x00000002u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000002u;
  }
  _impl_.version_.SetAllocated(version, GetArenaForAllocation());
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.async_fusion_config.ModelConfigs.version)
}

// optional string match_method = 3 [default = "hm_matcher"];
inline bool ModelConfigs::_internal_has_match_method() const {
  bool value = (_impl_._has_bits_[0] & 0x00000004u) != 0;
  return value;
}
inline bool ModelConfigs::has_match_method() const {
  return _internal_has_match_method();
}
inline void ModelConfigs::clear_match_method() {
  _impl_.match_method_.ClearToDefault(::apollo::perception::async_fusion_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_match_method_, GetArenaForAllocation());
  _impl_._has_bits_[0] &= ~0x00000004u;
}
inline const std::string& ModelConfigs::match_method() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.match_method)
  if (_impl_.match_method_.IsDefault()) return Impl_::_i_give_permission_to_break_this_code_default_match_method_.get();
  return _internal_match_method();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void ModelConfigs::set_match_method(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000004u;
 _impl_.match_method_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.match_method)
}
inline std::string* ModelConfigs::mutable_match_method() {
  std::string* _s = _internal_mutable_match_method();
  // @@protoc_insertion_point(field_mutable:apollo.perception.async_fusion_config.ModelConfigs.match_method)
  return _s;
}
inline const std::string& ModelConfigs::_internal_match_method() const {
  return _impl_.match_method_.Get();
}
inline void ModelConfigs::_internal_set_match_method(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000004u;
  _impl_.match_method_.Set(value, GetArenaForAllocation());
}
inline std::string* ModelConfigs::_internal_mutable_match_method() {
  _impl_._has_bits_[0] |= 0x00000004u;
  return _impl_.match_method_.Mutable(::apollo::perception::async_fusion_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_match_method_, GetArenaForAllocation());
}
inline std::string* ModelConfigs::release_match_method() {
  // @@protoc_insertion_point(field_release:apollo.perception.async_fusion_config.ModelConfigs.match_method)
  if (!_internal_has_match_method()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000004u;
  auto* p = _impl_.match_method_.Release();
  return p;
}
inline void ModelConfigs::set_allocated_match_method(std::string* match_method) {
  if (match_method != nullptr) {
    _impl_._has_bits_[0] |= 0x00000004u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000004u;
  }
  _impl_.match_method_.SetAllocated(match_method, GetArenaForAllocation());
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.async_fusion_config.ModelConfigs.match_method)
}

// optional float max_match_distance = 4 [default = 4];
inline bool ModelConfigs::_internal_has_max_match_distance() const {
  bool value = (_impl_._has_bits_[0] & 0x00000010u) != 0;
  return value;
}
inline bool ModelConfigs::has_max_match_distance() const {
  return _internal_has_max_match_distance();
}
inline void ModelConfigs::clear_max_match_distance() {
  _impl_.max_match_distance_ = 4;
  _impl_._has_bits_[0] &= ~0x00000010u;
}
inline float ModelConfigs::_internal_max_match_distance() const {
  return _impl_.max_match_distance_;
}
inline float ModelConfigs::max_match_distance() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.max_match_distance)
  return _internal_max_match_distance();
}
inline void ModelConfigs::_internal_set_max_match_distance(float value) {
  _impl_._has_bits_[0] |= 0x00000010u;
  _impl_.max_match_distance_ = value;
}
inline void ModelConfigs::set_max_match_distance(float value) {
  _internal_set_max_match_distance(value);
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.max_match_distance)
}

// optional float max_lidar_invisible_period = 5 [default = 0.25];
inline bool ModelConfigs::_internal_has_max_lidar_invisible_period() const {
  bool value = (_impl_._has_bits_[0] & 0x00000020u) != 0;
  return value;
}
inline bool ModelConfigs::has_max_lidar_invisible_period() const {
  return _internal_has_max_lidar_invisible_period();
}
inline void ModelConfigs::clear_max_lidar_invisible_period() {
  _impl_.max_lidar_invisible_period_ = 0.25f;
  _impl_._has_bits_[0] &= ~0x00000020u;
}
inline float ModelConfigs::_internal_max_lidar_invisible_period() const {
  return _impl_.max_lidar_invisible_period_;
}
inline float ModelConfigs::max_lidar_invisible_period() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.max_lidar_invisible_period)
  return _internal_max_lidar_invisible_period();
}
inline void ModelConfigs::_internal_set_max_lidar_invisible_period(float value) {
  _impl_._has_bits_[0] |= 0x00000020u;
  _impl_.max_lidar_invisible_period_ = value;
}
inline void ModelConfigs::set_max_lidar_invisible_period(float value) {
  _internal_set_max_lidar_invisible_period(value);
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.max_lidar_invisible_period)
}

// optional float max_radar_invisible_period = 6 [default = 0.25];
inline bool ModelConfigs::_internal_has_max_radar_invisible_period() const {
  bool value = (_impl_._has_bits_[0] & 0x00000040u) != 0;
  return value;
}
inline bool ModelConfigs::has_max_radar_invisible_period() const {
  return _internal_has_max_radar_invisible_period();
}
inline void ModelConfigs::clear_max_radar_invisible_period() {
  _impl_.max_radar_invisible_period_ = 0.25f;
  _impl_._has_bits_[0] &= ~0x00000040u;
}
inline float ModelConfigs::_internal_max_radar_invisible_period() const {
  return _impl_.max_radar_invisible_period_;
}
inline float ModelConfigs::max_radar_invisible_period() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.max_radar_invisible_period)
  return _internal_max_radar_invisible_period();
}
inline void ModelConfigs::_internal_set_max_radar_invisible_period(float value) {
  _impl_._has_bits_[0] |= 0x00000040u;
  _impl_.max_radar_invisible_period_ = value;
}
inline void ModelConfigs::set_max_radar_invisible_period(float value) {
  _internal_set_max_radar_invisible_period(value);
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.max_radar_invisible_period)
}

// optional float max_radar_confident_angle = 7 [default = 30];
inline bool ModelConfigs::_internal_has_max_radar_confident_angle() const {
  bool value = (_impl_._has_bits_[0] & 0x00000080u) != 0;
  return value;
}
inline bool ModelConfigs::has_max_radar_confident_angle() const {
  return _internal_has_max_radar_confident_angle();
}
inline void ModelConfigs::clear_max_radar_confident_angle() {
  _impl_.max_radar_confident_angle_ = 30;
  _impl_._has_bits_[0] &= ~0x00000080u;
}
inline float ModelConfigs::_internal_max_radar_confident_angle() const {
  return _impl_.max_radar_confident_angle_;
}
inline float ModelConfigs::max_radar_confident_angle() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.max_radar_confident_angle)
  return _internal_max_radar_confident_angle();
}
inline void ModelConfigs::_internal_set_max_radar_confident_angle(float value) {
  _impl_._has_bits_[0] |= 0x00000080u;
  _impl_.max_radar_confident_angle_ = value;
}
inline void ModelConfigs::set_max_radar_confident_angle(float value) {
  _internal_set_max_radar_confident_angle(value);
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.max_radar_confident_angle)
}

// optional float min_radar_confident_distance = 8 [default = 40];
inline bool ModelConfigs::_internal_has_min_radar_confident_distance() const {
  bool value = (_impl_._has_bits_[0] & 0x00000100u) != 0;
  return value;
}
inline bool ModelConfigs::has_min_radar_confident_distance() const {
  return _internal_has_min_radar_confident_distance();
}
inline void ModelConfigs::clear_min_radar_confident_distance() {
  _impl_.min_radar_confident_distance_ = 40;
  _impl_._has_bits_[0] &= ~0x00000100u;
}
inline float ModelConfigs::_internal_min_radar_confident_distance() const {
  return _impl_.min_radar_confident_distance_;
}
inline float ModelConfigs::min_radar_confident_distance() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.min_radar_confident_distance)
  return _internal_min_radar_confident_distance();
}
inline void ModelConfigs::_internal_set_min_radar_confident_distance(float value) {
  _impl_._has_bits_[0] |= 0x00000100u;
  _impl_.min_radar_confident_distance_ = value;
}
inline void ModelConfigs::set_min_radar_confident_distance(float value) {
  _internal_set_min_radar_confident_distance(value);
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.min_radar_confident_distance)
}

// optional bool publish_if_has_lidar = 9 [default = true];
inline bool ModelConfigs::_internal_has_publish_if_has_lidar() const {
  bool value = (_impl_._has_bits_[0] & 0x00000200u) != 0;
  return value;
}
inline bool ModelConfigs::has_publish_if_has_lidar() const {
  return _internal_has_publish_if_has_lidar();
}
inline void ModelConfigs::clear_publish_if_has_lidar() {
  _impl_.publish_if_has_lidar_ = true;
  _impl_._has_bits_[0] &= ~0x00000200u;
}
inline bool ModelConfigs::_internal_publish_if_has_lidar() const {
  return _impl_.publish_if_has_lidar_;
}
inline bool ModelConfigs::publish_if_has_lidar() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.publish_if_has_lidar)
  return _internal_publish_if_has_lidar();
}
inline void ModelConfigs::_internal_set_publish_if_has_lidar(bool value) {
  _impl_._has_bits_[0] |= 0x00000200u;
  _impl_.publish_if_has_lidar_ = value;
}
inline void ModelConfigs::set_publish_if_has_lidar(bool value) {
  _internal_set_publish_if_has_lidar(value);
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.publish_if_has_lidar)
}

// optional bool publish_if_has_radar = 10 [default = true];
inline bool ModelConfigs::_internal_has_publish_if_has_radar() const {
  bool value = (_impl_._has_bits_[0] & 0x00000400u) != 0;
  return value;
}
inline bool ModelConfigs::has_publish_if_has_radar() const {
  return _internal_has_publish_if_has_radar();
}
inline void ModelConfigs::clear_publish_if_has_radar() {
  _impl_.publish_if_has_radar_ = true;
  _impl_._has_bits_[0] &= ~0x00000400u;
}
inline bool ModelConfigs::_internal_publish_if_has_radar() const {
  return _impl_.publish_if_has_radar_;
}
inline bool ModelConfigs::publish_if_has_radar() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.publish_if_has_radar)
  return _internal_publish_if_has_radar();
}
inline void ModelConfigs::_internal_set_publish_if_has_radar(bool value) {
  _impl_._has_bits_[0] |= 0x00000400u;
  _impl_.publish_if_has_radar_ = value;
}
inline void ModelConfigs::set_publish_if_has_radar(bool value) {
  _internal_set_publish_if_has_radar(value);
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.publish_if_has_radar)
}

// optional string publish_sensor = 11 [default = "velodyne_64"];
inline bool ModelConfigs::_internal_has_publish_sensor() const {
  bool value = (_impl_._has_bits_[0] & 0x00000008u) != 0;
  return value;
}
inline bool ModelConfigs::has_publish_sensor() const {
  return _internal_has_publish_sensor();
}
inline void ModelConfigs::clear_publish_sensor() {
  _impl_.publish_sensor_.ClearToDefault(::apollo::perception::async_fusion_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_publish_sensor_, GetArenaForAllocation());
  _impl_._has_bits_[0] &= ~0x00000008u;
}
inline const std::string& ModelConfigs::publish_sensor() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.publish_sensor)
  if (_impl_.publish_sensor_.IsDefault()) return Impl_::_i_give_permission_to_break_this_code_default_publish_sensor_.get();
  return _internal_publish_sensor();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void ModelConfigs::set_publish_sensor(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000008u;
 _impl_.publish_sensor_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.publish_sensor)
}
inline std::string* ModelConfigs::mutable_publish_sensor() {
  std::string* _s = _internal_mutable_publish_sensor();
  // @@protoc_insertion_point(field_mutable:apollo.perception.async_fusion_config.ModelConfigs.publish_sensor)
  return _s;
}
inline const std::string& ModelConfigs::_internal_publish_sensor() const {
  return _impl_.publish_sensor_.Get();
}
inline void ModelConfigs::_internal_set_publish_sensor(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000008u;
  _impl_.publish_sensor_.Set(value, GetArenaForAllocation());
}
inline std::string* ModelConfigs::_internal_mutable_publish_sensor() {
  _impl_._has_bits_[0] |= 0x00000008u;
  return _impl_.publish_sensor_.Mutable(::apollo::perception::async_fusion_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_publish_sensor_, GetArenaForAllocation());
}
inline std::string* ModelConfigs::release_publish_sensor() {
  // @@protoc_insertion_point(field_release:apollo.perception.async_fusion_config.ModelConfigs.publish_sensor)
  if (!_internal_has_publish_sensor()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000008u;
  auto* p = _impl_.publish_sensor_.Release();
  return p;
}
inline void ModelConfigs::set_allocated_publish_sensor(std::string* publish_sensor) {
  if (publish_sensor != nullptr) {
    _impl_._has_bits_[0] |= 0x00000008u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000008u;
  }
  _impl_.publish_sensor_.SetAllocated(publish_sensor, GetArenaForAllocation());
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.async_fusion_config.ModelConfigs.publish_sensor)
}

// optional bool use_radar = 12 [default = true];
inline bool ModelConfigs::_internal_has_use_radar() const {
  bool value = (_impl_._has_bits_[0] & 0x00000800u) != 0;
  return value;
}
inline bool ModelConfigs::has_use_radar() const {
  return _internal_has_use_radar();
}
inline void ModelConfigs::clear_use_radar() {
  _impl_.use_radar_ = true;
  _impl_._has_bits_[0] &= ~0x00000800u;
}
inline bool ModelConfigs::_internal_use_radar() const {
  return _impl_.use_radar_;
}
inline bool ModelConfigs::use_radar() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.use_radar)
  return _internal_use_radar();
}
inline void ModelConfigs::_internal_set_use_radar(bool value) {
  _impl_._has_bits_[0] |= 0x00000800u;
  _impl_.use_radar_ = value;
}
inline void ModelConfigs::set_use_radar(bool value) {
  _internal_set_use_radar(value);
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.use_radar)
}

// optional bool use_lidar = 13 [default = true];
inline bool ModelConfigs::_internal_has_use_lidar() const {
  bool value = (_impl_._has_bits_[0] & 0x00001000u) != 0;
  return value;
}
inline bool ModelConfigs::has_use_lidar() const {
  return _internal_has_use_lidar();
}
inline void ModelConfigs::clear_use_lidar() {
  _impl_.use_lidar_ = true;
  _impl_._has_bits_[0] &= ~0x00001000u;
}
inline bool ModelConfigs::_internal_use_lidar() const {
  return _impl_.use_lidar_;
}
inline bool ModelConfigs::use_lidar() const {
  // @@protoc_insertion_point(field_get:apollo.perception.async_fusion_config.ModelConfigs.use_lidar)
  return _internal_use_lidar();
}
inline void ModelConfigs::_internal_set_use_lidar(bool value) {
  _impl_._has_bits_[0] |= 0x00001000u;
  _impl_.use_lidar_ = value;
}
inline void ModelConfigs::set_use_lidar(bool value) {
  _internal_set_use_lidar(value);
  // @@protoc_insertion_point(field_set:apollo.perception.async_fusion_config.ModelConfigs.use_lidar)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace async_fusion_config
}  // namespace perception
}  // namespace apollo

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fproto_2fasync_5ffusion_5fconfig_2eproto
