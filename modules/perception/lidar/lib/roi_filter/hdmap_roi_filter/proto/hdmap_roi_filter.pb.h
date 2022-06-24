// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/perception/lidar/lib/roi_filter/hdmap_roi_filter/proto/hdmap_roi_filter.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2flidar_2flib_2froi_5ffilter_2fhdmap_5froi_5ffilter_2fproto_2fhdmap_5froi_5ffilter_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2flidar_2flib_2froi_5ffilter_2fhdmap_5froi_5ffilter_2fproto_2fhdmap_5froi_5ffilter_2eproto

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
#define PROTOBUF_INTERNAL_EXPORT_modules_2fperception_2flidar_2flib_2froi_5ffilter_2fhdmap_5froi_5ffilter_2fproto_2fhdmap_5froi_5ffilter_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_modules_2fperception_2flidar_2flib_2froi_5ffilter_2fhdmap_5froi_5ffilter_2fproto_2fhdmap_5froi_5ffilter_2eproto {
  static const uint32_t offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_modules_2fperception_2flidar_2flib_2froi_5ffilter_2fhdmap_5froi_5ffilter_2fproto_2fhdmap_5froi_5ffilter_2eproto;
namespace apollo {
namespace perception {
namespace lidar {
class HDMapRoiFilterConfig;
struct HDMapRoiFilterConfigDefaultTypeInternal;
extern HDMapRoiFilterConfigDefaultTypeInternal _HDMapRoiFilterConfig_default_instance_;
}  // namespace lidar
}  // namespace perception
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> ::apollo::perception::lidar::HDMapRoiFilterConfig* Arena::CreateMaybeMessage<::apollo::perception::lidar::HDMapRoiFilterConfig>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace apollo {
namespace perception {
namespace lidar {

// ===================================================================

class HDMapRoiFilterConfig final :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:apollo.perception.lidar.HDMapRoiFilterConfig) */ {
 public:
  inline HDMapRoiFilterConfig() : HDMapRoiFilterConfig(nullptr) {}
  ~HDMapRoiFilterConfig() override;
  explicit PROTOBUF_CONSTEXPR HDMapRoiFilterConfig(::PROTOBUF_NAMESPACE_ID::internal::ConstantInitialized);

  HDMapRoiFilterConfig(const HDMapRoiFilterConfig& from);
  HDMapRoiFilterConfig(HDMapRoiFilterConfig&& from) noexcept
    : HDMapRoiFilterConfig() {
    *this = ::std::move(from);
  }

  inline HDMapRoiFilterConfig& operator=(const HDMapRoiFilterConfig& from) {
    CopyFrom(from);
    return *this;
  }
  inline HDMapRoiFilterConfig& operator=(HDMapRoiFilterConfig&& from) noexcept {
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
  static const HDMapRoiFilterConfig& default_instance() {
    return *internal_default_instance();
  }
  static inline const HDMapRoiFilterConfig* internal_default_instance() {
    return reinterpret_cast<const HDMapRoiFilterConfig*>(
               &_HDMapRoiFilterConfig_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(HDMapRoiFilterConfig& a, HDMapRoiFilterConfig& b) {
    a.Swap(&b);
  }
  inline void Swap(HDMapRoiFilterConfig* other) {
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
  void UnsafeArenaSwap(HDMapRoiFilterConfig* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetOwningArena() == other->GetOwningArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  HDMapRoiFilterConfig* New(::PROTOBUF_NAMESPACE_ID::Arena* arena = nullptr) const final {
    return CreateMaybeMessage<HDMapRoiFilterConfig>(arena);
  }
  using ::PROTOBUF_NAMESPACE_ID::Message::CopyFrom;
  void CopyFrom(const HDMapRoiFilterConfig& from);
  using ::PROTOBUF_NAMESPACE_ID::Message::MergeFrom;
  void MergeFrom(const HDMapRoiFilterConfig& from);
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
  void InternalSwap(HDMapRoiFilterConfig* other);

  private:
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "apollo.perception.lidar.HDMapRoiFilterConfig";
  }
  protected:
  explicit HDMapRoiFilterConfig(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                       bool is_message_owned = false);
  public:

  static const ClassData _class_data_;
  const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*GetClassData() const final;

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kExtendDistFieldNumber = 3,
    kNoEdgeTableFieldNumber = 4,
    kSetRoiServiceFieldNumber = 5,
    kRangeFieldNumber = 1,
    kCellSizeFieldNumber = 2,
  };
  // optional double extend_dist = 3 [default = 0];
  bool has_extend_dist() const;
  private:
  bool _internal_has_extend_dist() const;
  public:
  void clear_extend_dist();
  double extend_dist() const;
  void set_extend_dist(double value);
  private:
  double _internal_extend_dist() const;
  void _internal_set_extend_dist(double value);
  public:

  // optional bool no_edge_table = 4 [default = false];
  bool has_no_edge_table() const;
  private:
  bool _internal_has_no_edge_table() const;
  public:
  void clear_no_edge_table();
  bool no_edge_table() const;
  void set_no_edge_table(bool value);
  private:
  bool _internal_no_edge_table() const;
  void _internal_set_no_edge_table(bool value);
  public:

  // optional bool set_roi_service = 5 [default = false];
  bool has_set_roi_service() const;
  private:
  bool _internal_has_set_roi_service() const;
  public:
  void clear_set_roi_service();
  bool set_roi_service() const;
  void set_set_roi_service(bool value);
  private:
  bool _internal_set_roi_service() const;
  void _internal_set_set_roi_service(bool value);
  public:

  // optional double range = 1 [default = 120];
  bool has_range() const;
  private:
  bool _internal_has_range() const;
  public:
  void clear_range();
  double range() const;
  void set_range(double value);
  private:
  double _internal_range() const;
  void _internal_set_range(double value);
  public:

  // optional double cell_size = 2 [default = 0.25];
  bool has_cell_size() const;
  private:
  bool _internal_has_cell_size() const;
  public:
  void clear_cell_size();
  double cell_size() const;
  void set_cell_size(double value);
  private:
  double _internal_cell_size() const;
  void _internal_set_cell_size(double value);
  public:

  // @@protoc_insertion_point(class_scope:apollo.perception.lidar.HDMapRoiFilterConfig)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  struct Impl_ {
    ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
    mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
    double extend_dist_;
    bool no_edge_table_;
    bool set_roi_service_;
    double range_;
    double cell_size_;
  };
  union { Impl_ _impl_; };
  friend struct ::TableStruct_modules_2fperception_2flidar_2flib_2froi_5ffilter_2fhdmap_5froi_5ffilter_2fproto_2fhdmap_5froi_5ffilter_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// HDMapRoiFilterConfig

// optional double range = 1 [default = 120];
inline bool HDMapRoiFilterConfig::_internal_has_range() const {
  bool value = (_impl_._has_bits_[0] & 0x00000008u) != 0;
  return value;
}
inline bool HDMapRoiFilterConfig::has_range() const {
  return _internal_has_range();
}
inline void HDMapRoiFilterConfig::clear_range() {
  _impl_.range_ = 120;
  _impl_._has_bits_[0] &= ~0x00000008u;
}
inline double HDMapRoiFilterConfig::_internal_range() const {
  return _impl_.range_;
}
inline double HDMapRoiFilterConfig::range() const {
  // @@protoc_insertion_point(field_get:apollo.perception.lidar.HDMapRoiFilterConfig.range)
  return _internal_range();
}
inline void HDMapRoiFilterConfig::_internal_set_range(double value) {
  _impl_._has_bits_[0] |= 0x00000008u;
  _impl_.range_ = value;
}
inline void HDMapRoiFilterConfig::set_range(double value) {
  _internal_set_range(value);
  // @@protoc_insertion_point(field_set:apollo.perception.lidar.HDMapRoiFilterConfig.range)
}

// optional double cell_size = 2 [default = 0.25];
inline bool HDMapRoiFilterConfig::_internal_has_cell_size() const {
  bool value = (_impl_._has_bits_[0] & 0x00000010u) != 0;
  return value;
}
inline bool HDMapRoiFilterConfig::has_cell_size() const {
  return _internal_has_cell_size();
}
inline void HDMapRoiFilterConfig::clear_cell_size() {
  _impl_.cell_size_ = 0.25;
  _impl_._has_bits_[0] &= ~0x00000010u;
}
inline double HDMapRoiFilterConfig::_internal_cell_size() const {
  return _impl_.cell_size_;
}
inline double HDMapRoiFilterConfig::cell_size() const {
  // @@protoc_insertion_point(field_get:apollo.perception.lidar.HDMapRoiFilterConfig.cell_size)
  return _internal_cell_size();
}
inline void HDMapRoiFilterConfig::_internal_set_cell_size(double value) {
  _impl_._has_bits_[0] |= 0x00000010u;
  _impl_.cell_size_ = value;
}
inline void HDMapRoiFilterConfig::set_cell_size(double value) {
  _internal_set_cell_size(value);
  // @@protoc_insertion_point(field_set:apollo.perception.lidar.HDMapRoiFilterConfig.cell_size)
}

// optional double extend_dist = 3 [default = 0];
inline bool HDMapRoiFilterConfig::_internal_has_extend_dist() const {
  bool value = (_impl_._has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool HDMapRoiFilterConfig::has_extend_dist() const {
  return _internal_has_extend_dist();
}
inline void HDMapRoiFilterConfig::clear_extend_dist() {
  _impl_.extend_dist_ = 0;
  _impl_._has_bits_[0] &= ~0x00000001u;
}
inline double HDMapRoiFilterConfig::_internal_extend_dist() const {
  return _impl_.extend_dist_;
}
inline double HDMapRoiFilterConfig::extend_dist() const {
  // @@protoc_insertion_point(field_get:apollo.perception.lidar.HDMapRoiFilterConfig.extend_dist)
  return _internal_extend_dist();
}
inline void HDMapRoiFilterConfig::_internal_set_extend_dist(double value) {
  _impl_._has_bits_[0] |= 0x00000001u;
  _impl_.extend_dist_ = value;
}
inline void HDMapRoiFilterConfig::set_extend_dist(double value) {
  _internal_set_extend_dist(value);
  // @@protoc_insertion_point(field_set:apollo.perception.lidar.HDMapRoiFilterConfig.extend_dist)
}

// optional bool no_edge_table = 4 [default = false];
inline bool HDMapRoiFilterConfig::_internal_has_no_edge_table() const {
  bool value = (_impl_._has_bits_[0] & 0x00000002u) != 0;
  return value;
}
inline bool HDMapRoiFilterConfig::has_no_edge_table() const {
  return _internal_has_no_edge_table();
}
inline void HDMapRoiFilterConfig::clear_no_edge_table() {
  _impl_.no_edge_table_ = false;
  _impl_._has_bits_[0] &= ~0x00000002u;
}
inline bool HDMapRoiFilterConfig::_internal_no_edge_table() const {
  return _impl_.no_edge_table_;
}
inline bool HDMapRoiFilterConfig::no_edge_table() const {
  // @@protoc_insertion_point(field_get:apollo.perception.lidar.HDMapRoiFilterConfig.no_edge_table)
  return _internal_no_edge_table();
}
inline void HDMapRoiFilterConfig::_internal_set_no_edge_table(bool value) {
  _impl_._has_bits_[0] |= 0x00000002u;
  _impl_.no_edge_table_ = value;
}
inline void HDMapRoiFilterConfig::set_no_edge_table(bool value) {
  _internal_set_no_edge_table(value);
  // @@protoc_insertion_point(field_set:apollo.perception.lidar.HDMapRoiFilterConfig.no_edge_table)
}

// optional bool set_roi_service = 5 [default = false];
inline bool HDMapRoiFilterConfig::_internal_has_set_roi_service() const {
  bool value = (_impl_._has_bits_[0] & 0x00000004u) != 0;
  return value;
}
inline bool HDMapRoiFilterConfig::has_set_roi_service() const {
  return _internal_has_set_roi_service();
}
inline void HDMapRoiFilterConfig::clear_set_roi_service() {
  _impl_.set_roi_service_ = false;
  _impl_._has_bits_[0] &= ~0x00000004u;
}
inline bool HDMapRoiFilterConfig::_internal_set_roi_service() const {
  return _impl_.set_roi_service_;
}
inline bool HDMapRoiFilterConfig::set_roi_service() const {
  // @@protoc_insertion_point(field_get:apollo.perception.lidar.HDMapRoiFilterConfig.set_roi_service)
  return _internal_set_roi_service();
}
inline void HDMapRoiFilterConfig::_internal_set_set_roi_service(bool value) {
  _impl_._has_bits_[0] |= 0x00000004u;
  _impl_.set_roi_service_ = value;
}
inline void HDMapRoiFilterConfig::set_set_roi_service(bool value) {
  _internal_set_set_roi_service(value);
  // @@protoc_insertion_point(field_set:apollo.perception.lidar.HDMapRoiFilterConfig.set_roi_service)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace lidar
}  // namespace perception
}  // namespace apollo

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2flidar_2flib_2froi_5ffilter_2fhdmap_5froi_5ffilter_2fproto_2fhdmap_5froi_5ffilter_2eproto
