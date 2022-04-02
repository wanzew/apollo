// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/perception/proto/cnn_segmentation_config.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fproto_2fcnn_5fsegmentation_5fconfig_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fproto_2fcnn_5fsegmentation_5fconfig_2eproto

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
#define PROTOBUF_INTERNAL_EXPORT_modules_2fperception_2fproto_2fcnn_5fsegmentation_5fconfig_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_modules_2fperception_2fproto_2fcnn_5fsegmentation_5fconfig_2eproto {
  static const uint32_t offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_modules_2fperception_2fproto_2fcnn_5fsegmentation_5fconfig_2eproto;
namespace apollo {
namespace perception {
namespace cnn_segmentation_config {
class ModelConfigs;
struct ModelConfigsDefaultTypeInternal;
extern ModelConfigsDefaultTypeInternal _ModelConfigs_default_instance_;
}  // namespace cnn_segmentation_config
}  // namespace perception
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> ::apollo::perception::cnn_segmentation_config::ModelConfigs* Arena::CreateMaybeMessage<::apollo::perception::cnn_segmentation_config::ModelConfigs>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace apollo {
namespace perception {
namespace cnn_segmentation_config {

// ===================================================================

class ModelConfigs final :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:apollo.perception.cnn_segmentation_config.ModelConfigs) */ {
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
    return "apollo.perception.cnn_segmentation_config.ModelConfigs";
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
    kConfigFileFieldNumber = 3,
    kProtoFileFieldNumber = 4,
    kWeightFileFieldNumber = 5,
  };
  // optional string name = 1 [default = "CNNSegmentation"];
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

  // optional string config_file = 3 [default = "modules/perception/model/cnn_segmentation/cnnseg.conf"];
  bool has_config_file() const;
  private:
  bool _internal_has_config_file() const;
  public:
  void clear_config_file();
  const std::string& config_file() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_config_file(ArgT0&& arg0, ArgT... args);
  std::string* mutable_config_file();
  PROTOBUF_NODISCARD std::string* release_config_file();
  void set_allocated_config_file(std::string* config_file);
  private:
  const std::string& _internal_config_file() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_config_file(const std::string& value);
  std::string* _internal_mutable_config_file();
  public:

  // optional string proto_file = 4 [default = "modules/perception/model/cnn_segmentation/deploy.prototxt"];
  bool has_proto_file() const;
  private:
  bool _internal_has_proto_file() const;
  public:
  void clear_proto_file();
  const std::string& proto_file() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_proto_file(ArgT0&& arg0, ArgT... args);
  std::string* mutable_proto_file();
  PROTOBUF_NODISCARD std::string* release_proto_file();
  void set_allocated_proto_file(std::string* proto_file);
  private:
  const std::string& _internal_proto_file() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_proto_file(const std::string& value);
  std::string* _internal_mutable_proto_file();
  public:

  // optional string weight_file = 5 [default = "modules/perception/model/cnn_segmentation/deploy.caffemodel"];
  bool has_weight_file() const;
  private:
  bool _internal_has_weight_file() const;
  public:
  void clear_weight_file();
  const std::string& weight_file() const;
  template <typename ArgT0 = const std::string&, typename... ArgT>
  void set_weight_file(ArgT0&& arg0, ArgT... args);
  std::string* mutable_weight_file();
  PROTOBUF_NODISCARD std::string* release_weight_file();
  void set_allocated_weight_file(std::string* weight_file);
  private:
  const std::string& _internal_weight_file() const;
  inline PROTOBUF_ALWAYS_INLINE void _internal_set_weight_file(const std::string& value);
  std::string* _internal_mutable_weight_file();
  public:

  // @@protoc_insertion_point(class_scope:apollo.perception.cnn_segmentation_config.ModelConfigs)
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
    static const ::PROTOBUF_NAMESPACE_ID::internal::LazyString _i_give_permission_to_break_this_code_default_config_file_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr config_file_;
    static const ::PROTOBUF_NAMESPACE_ID::internal::LazyString _i_give_permission_to_break_this_code_default_proto_file_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr proto_file_;
    static const ::PROTOBUF_NAMESPACE_ID::internal::LazyString _i_give_permission_to_break_this_code_default_weight_file_;
    ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr weight_file_;
  };
  union { Impl_ _impl_; };
  friend struct ::TableStruct_modules_2fperception_2fproto_2fcnn_5fsegmentation_5fconfig_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// ModelConfigs

// optional string name = 1 [default = "CNNSegmentation"];
inline bool ModelConfigs::_internal_has_name() const {
  bool value = (_impl_._has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool ModelConfigs::has_name() const {
  return _internal_has_name();
}
inline void ModelConfigs::clear_name() {
  _impl_.name_.ClearToDefault(::apollo::perception::cnn_segmentation_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_name_, GetArenaForAllocation());
  _impl_._has_bits_[0] &= ~0x00000001u;
}
inline const std::string& ModelConfigs::name() const {
  // @@protoc_insertion_point(field_get:apollo.perception.cnn_segmentation_config.ModelConfigs.name)
  if (_impl_.name_.IsDefault()) return Impl_::_i_give_permission_to_break_this_code_default_name_.get();
  return _internal_name();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void ModelConfigs::set_name(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000001u;
 _impl_.name_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.cnn_segmentation_config.ModelConfigs.name)
}
inline std::string* ModelConfigs::mutable_name() {
  std::string* _s = _internal_mutable_name();
  // @@protoc_insertion_point(field_mutable:apollo.perception.cnn_segmentation_config.ModelConfigs.name)
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
  return _impl_.name_.Mutable(::apollo::perception::cnn_segmentation_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_name_, GetArenaForAllocation());
}
inline std::string* ModelConfigs::release_name() {
  // @@protoc_insertion_point(field_release:apollo.perception.cnn_segmentation_config.ModelConfigs.name)
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
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.cnn_segmentation_config.ModelConfigs.name)
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
  _impl_.version_.ClearToDefault(::apollo::perception::cnn_segmentation_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_version_, GetArenaForAllocation());
  _impl_._has_bits_[0] &= ~0x00000002u;
}
inline const std::string& ModelConfigs::version() const {
  // @@protoc_insertion_point(field_get:apollo.perception.cnn_segmentation_config.ModelConfigs.version)
  if (_impl_.version_.IsDefault()) return Impl_::_i_give_permission_to_break_this_code_default_version_.get();
  return _internal_version();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void ModelConfigs::set_version(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000002u;
 _impl_.version_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.cnn_segmentation_config.ModelConfigs.version)
}
inline std::string* ModelConfigs::mutable_version() {
  std::string* _s = _internal_mutable_version();
  // @@protoc_insertion_point(field_mutable:apollo.perception.cnn_segmentation_config.ModelConfigs.version)
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
  return _impl_.version_.Mutable(::apollo::perception::cnn_segmentation_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_version_, GetArenaForAllocation());
}
inline std::string* ModelConfigs::release_version() {
  // @@protoc_insertion_point(field_release:apollo.perception.cnn_segmentation_config.ModelConfigs.version)
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
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.cnn_segmentation_config.ModelConfigs.version)
}

// optional string config_file = 3 [default = "modules/perception/model/cnn_segmentation/cnnseg.conf"];
inline bool ModelConfigs::_internal_has_config_file() const {
  bool value = (_impl_._has_bits_[0] & 0x00000004u) != 0;
  return value;
}
inline bool ModelConfigs::has_config_file() const {
  return _internal_has_config_file();
}
inline void ModelConfigs::clear_config_file() {
  _impl_.config_file_.ClearToDefault(::apollo::perception::cnn_segmentation_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_config_file_, GetArenaForAllocation());
  _impl_._has_bits_[0] &= ~0x00000004u;
}
inline const std::string& ModelConfigs::config_file() const {
  // @@protoc_insertion_point(field_get:apollo.perception.cnn_segmentation_config.ModelConfigs.config_file)
  if (_impl_.config_file_.IsDefault()) return Impl_::_i_give_permission_to_break_this_code_default_config_file_.get();
  return _internal_config_file();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void ModelConfigs::set_config_file(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000004u;
 _impl_.config_file_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.cnn_segmentation_config.ModelConfigs.config_file)
}
inline std::string* ModelConfigs::mutable_config_file() {
  std::string* _s = _internal_mutable_config_file();
  // @@protoc_insertion_point(field_mutable:apollo.perception.cnn_segmentation_config.ModelConfigs.config_file)
  return _s;
}
inline const std::string& ModelConfigs::_internal_config_file() const {
  return _impl_.config_file_.Get();
}
inline void ModelConfigs::_internal_set_config_file(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000004u;
  _impl_.config_file_.Set(value, GetArenaForAllocation());
}
inline std::string* ModelConfigs::_internal_mutable_config_file() {
  _impl_._has_bits_[0] |= 0x00000004u;
  return _impl_.config_file_.Mutable(::apollo::perception::cnn_segmentation_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_config_file_, GetArenaForAllocation());
}
inline std::string* ModelConfigs::release_config_file() {
  // @@protoc_insertion_point(field_release:apollo.perception.cnn_segmentation_config.ModelConfigs.config_file)
  if (!_internal_has_config_file()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000004u;
  auto* p = _impl_.config_file_.Release();
  return p;
}
inline void ModelConfigs::set_allocated_config_file(std::string* config_file) {
  if (config_file != nullptr) {
    _impl_._has_bits_[0] |= 0x00000004u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000004u;
  }
  _impl_.config_file_.SetAllocated(config_file, GetArenaForAllocation());
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.cnn_segmentation_config.ModelConfigs.config_file)
}

// optional string proto_file = 4 [default = "modules/perception/model/cnn_segmentation/deploy.prototxt"];
inline bool ModelConfigs::_internal_has_proto_file() const {
  bool value = (_impl_._has_bits_[0] & 0x00000008u) != 0;
  return value;
}
inline bool ModelConfigs::has_proto_file() const {
  return _internal_has_proto_file();
}
inline void ModelConfigs::clear_proto_file() {
  _impl_.proto_file_.ClearToDefault(::apollo::perception::cnn_segmentation_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_proto_file_, GetArenaForAllocation());
  _impl_._has_bits_[0] &= ~0x00000008u;
}
inline const std::string& ModelConfigs::proto_file() const {
  // @@protoc_insertion_point(field_get:apollo.perception.cnn_segmentation_config.ModelConfigs.proto_file)
  if (_impl_.proto_file_.IsDefault()) return Impl_::_i_give_permission_to_break_this_code_default_proto_file_.get();
  return _internal_proto_file();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void ModelConfigs::set_proto_file(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000008u;
 _impl_.proto_file_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.cnn_segmentation_config.ModelConfigs.proto_file)
}
inline std::string* ModelConfigs::mutable_proto_file() {
  std::string* _s = _internal_mutable_proto_file();
  // @@protoc_insertion_point(field_mutable:apollo.perception.cnn_segmentation_config.ModelConfigs.proto_file)
  return _s;
}
inline const std::string& ModelConfigs::_internal_proto_file() const {
  return _impl_.proto_file_.Get();
}
inline void ModelConfigs::_internal_set_proto_file(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000008u;
  _impl_.proto_file_.Set(value, GetArenaForAllocation());
}
inline std::string* ModelConfigs::_internal_mutable_proto_file() {
  _impl_._has_bits_[0] |= 0x00000008u;
  return _impl_.proto_file_.Mutable(::apollo::perception::cnn_segmentation_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_proto_file_, GetArenaForAllocation());
}
inline std::string* ModelConfigs::release_proto_file() {
  // @@protoc_insertion_point(field_release:apollo.perception.cnn_segmentation_config.ModelConfigs.proto_file)
  if (!_internal_has_proto_file()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000008u;
  auto* p = _impl_.proto_file_.Release();
  return p;
}
inline void ModelConfigs::set_allocated_proto_file(std::string* proto_file) {
  if (proto_file != nullptr) {
    _impl_._has_bits_[0] |= 0x00000008u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000008u;
  }
  _impl_.proto_file_.SetAllocated(proto_file, GetArenaForAllocation());
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.cnn_segmentation_config.ModelConfigs.proto_file)
}

// optional string weight_file = 5 [default = "modules/perception/model/cnn_segmentation/deploy.caffemodel"];
inline bool ModelConfigs::_internal_has_weight_file() const {
  bool value = (_impl_._has_bits_[0] & 0x00000010u) != 0;
  return value;
}
inline bool ModelConfigs::has_weight_file() const {
  return _internal_has_weight_file();
}
inline void ModelConfigs::clear_weight_file() {
  _impl_.weight_file_.ClearToDefault(::apollo::perception::cnn_segmentation_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_weight_file_, GetArenaForAllocation());
  _impl_._has_bits_[0] &= ~0x00000010u;
}
inline const std::string& ModelConfigs::weight_file() const {
  // @@protoc_insertion_point(field_get:apollo.perception.cnn_segmentation_config.ModelConfigs.weight_file)
  if (_impl_.weight_file_.IsDefault()) return Impl_::_i_give_permission_to_break_this_code_default_weight_file_.get();
  return _internal_weight_file();
}
template <typename ArgT0, typename... ArgT>
inline PROTOBUF_ALWAYS_INLINE
void ModelConfigs::set_weight_file(ArgT0&& arg0, ArgT... args) {
 _impl_._has_bits_[0] |= 0x00000010u;
 _impl_.weight_file_.Set(static_cast<ArgT0 &&>(arg0), args..., GetArenaForAllocation());
  // @@protoc_insertion_point(field_set:apollo.perception.cnn_segmentation_config.ModelConfigs.weight_file)
}
inline std::string* ModelConfigs::mutable_weight_file() {
  std::string* _s = _internal_mutable_weight_file();
  // @@protoc_insertion_point(field_mutable:apollo.perception.cnn_segmentation_config.ModelConfigs.weight_file)
  return _s;
}
inline const std::string& ModelConfigs::_internal_weight_file() const {
  return _impl_.weight_file_.Get();
}
inline void ModelConfigs::_internal_set_weight_file(const std::string& value) {
  _impl_._has_bits_[0] |= 0x00000010u;
  _impl_.weight_file_.Set(value, GetArenaForAllocation());
}
inline std::string* ModelConfigs::_internal_mutable_weight_file() {
  _impl_._has_bits_[0] |= 0x00000010u;
  return _impl_.weight_file_.Mutable(::apollo::perception::cnn_segmentation_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_weight_file_, GetArenaForAllocation());
}
inline std::string* ModelConfigs::release_weight_file() {
  // @@protoc_insertion_point(field_release:apollo.perception.cnn_segmentation_config.ModelConfigs.weight_file)
  if (!_internal_has_weight_file()) {
    return nullptr;
  }
  _impl_._has_bits_[0] &= ~0x00000010u;
  auto* p = _impl_.weight_file_.Release();
  return p;
}
inline void ModelConfigs::set_allocated_weight_file(std::string* weight_file) {
  if (weight_file != nullptr) {
    _impl_._has_bits_[0] |= 0x00000010u;
  } else {
    _impl_._has_bits_[0] &= ~0x00000010u;
  }
  _impl_.weight_file_.SetAllocated(weight_file, GetArenaForAllocation());
  // @@protoc_insertion_point(field_set_allocated:apollo.perception.cnn_segmentation_config.ModelConfigs.weight_file)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace cnn_segmentation_config
}  // namespace perception
}  // namespace apollo

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fproto_2fcnn_5fsegmentation_5fconfig_2eproto
