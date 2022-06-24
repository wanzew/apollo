// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/contrib/lgsvl_msgs/proto/detection3darray.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_modules_2fcontrib_2flgsvl_5fmsgs_2fproto_2fdetection3darray_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_modules_2fcontrib_2flgsvl_5fmsgs_2fproto_2fdetection3darray_2eproto

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
#include "modules/common/proto/header.pb.h"
#include "modules/contrib/lgsvl_msgs/proto/detection3d.pb.h"
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_modules_2fcontrib_2flgsvl_5fmsgs_2fproto_2fdetection3darray_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_modules_2fcontrib_2flgsvl_5fmsgs_2fproto_2fdetection3darray_2eproto {
  static const uint32_t offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_modules_2fcontrib_2flgsvl_5fmsgs_2fproto_2fdetection3darray_2eproto;
namespace apollo {
namespace contrib {
namespace lgsvl_msgs {
class Detection3DArray;
struct Detection3DArrayDefaultTypeInternal;
extern Detection3DArrayDefaultTypeInternal _Detection3DArray_default_instance_;
}  // namespace lgsvl_msgs
}  // namespace contrib
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> ::apollo::contrib::lgsvl_msgs::Detection3DArray* Arena::CreateMaybeMessage<::apollo::contrib::lgsvl_msgs::Detection3DArray>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace apollo {
namespace contrib {
namespace lgsvl_msgs {

// ===================================================================

class Detection3DArray final :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:apollo.contrib.lgsvl_msgs.Detection3DArray) */ {
 public:
  inline Detection3DArray() : Detection3DArray(nullptr) {}
  ~Detection3DArray() override;
  explicit PROTOBUF_CONSTEXPR Detection3DArray(::PROTOBUF_NAMESPACE_ID::internal::ConstantInitialized);

  Detection3DArray(const Detection3DArray& from);
  Detection3DArray(Detection3DArray&& from) noexcept
    : Detection3DArray() {
    *this = ::std::move(from);
  }

  inline Detection3DArray& operator=(const Detection3DArray& from) {
    CopyFrom(from);
    return *this;
  }
  inline Detection3DArray& operator=(Detection3DArray&& from) noexcept {
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

  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* descriptor() {
    return GetDescriptor();
  }
  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* GetDescriptor() {
    return default_instance().GetMetadata().descriptor;
  }
  static const ::PROTOBUF_NAMESPACE_ID::Reflection* GetReflection() {
    return default_instance().GetMetadata().reflection;
  }
  static const Detection3DArray& default_instance() {
    return *internal_default_instance();
  }
  static inline const Detection3DArray* internal_default_instance() {
    return reinterpret_cast<const Detection3DArray*>(
               &_Detection3DArray_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(Detection3DArray& a, Detection3DArray& b) {
    a.Swap(&b);
  }
  inline void Swap(Detection3DArray* other) {
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
  void UnsafeArenaSwap(Detection3DArray* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetOwningArena() == other->GetOwningArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  Detection3DArray* New(::PROTOBUF_NAMESPACE_ID::Arena* arena = nullptr) const final {
    return CreateMaybeMessage<Detection3DArray>(arena);
  }
  using ::PROTOBUF_NAMESPACE_ID::Message::CopyFrom;
  void CopyFrom(const Detection3DArray& from);
  using ::PROTOBUF_NAMESPACE_ID::Message::MergeFrom;
  void MergeFrom(const Detection3DArray& from);
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
  void InternalSwap(Detection3DArray* other);

  private:
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "apollo.contrib.lgsvl_msgs.Detection3DArray";
  }
  protected:
  explicit Detection3DArray(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                       bool is_message_owned = false);
  public:

  static const ClassData _class_data_;
  const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*GetClassData() const final;

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kDetectionsFieldNumber = 2,
    kHeaderFieldNumber = 1,
  };
  // repeated .apollo.contrib.lgsvl_msgs.Detection3D detections = 2;
  int detections_size() const;
  private:
  int _internal_detections_size() const;
  public:
  void clear_detections();
  ::apollo::contrib::lgsvl_msgs::Detection3D* mutable_detections(int index);
  ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::apollo::contrib::lgsvl_msgs::Detection3D >*
      mutable_detections();
  private:
  const ::apollo::contrib::lgsvl_msgs::Detection3D& _internal_detections(int index) const;
  ::apollo::contrib::lgsvl_msgs::Detection3D* _internal_add_detections();
  public:
  const ::apollo::contrib::lgsvl_msgs::Detection3D& detections(int index) const;
  ::apollo::contrib::lgsvl_msgs::Detection3D* add_detections();
  const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::apollo::contrib::lgsvl_msgs::Detection3D >&
      detections() const;

  // .apollo.common.Header header = 1;
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

  // @@protoc_insertion_point(class_scope:apollo.contrib.lgsvl_msgs.Detection3DArray)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  struct Impl_ {
    ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::apollo::contrib::lgsvl_msgs::Detection3D > detections_;
    ::apollo::common::Header* header_;
    mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  };
  union { Impl_ _impl_; };
  friend struct ::TableStruct_modules_2fcontrib_2flgsvl_5fmsgs_2fproto_2fdetection3darray_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// Detection3DArray

// .apollo.common.Header header = 1;
inline bool Detection3DArray::_internal_has_header() const {
  return this != internal_default_instance() && _impl_.header_ != nullptr;
}
inline bool Detection3DArray::has_header() const {
  return _internal_has_header();
}
inline const ::apollo::common::Header& Detection3DArray::_internal_header() const {
  const ::apollo::common::Header* p = _impl_.header_;
  return p != nullptr ? *p : reinterpret_cast<const ::apollo::common::Header&>(
      ::apollo::common::_Header_default_instance_);
}
inline const ::apollo::common::Header& Detection3DArray::header() const {
  // @@protoc_insertion_point(field_get:apollo.contrib.lgsvl_msgs.Detection3DArray.header)
  return _internal_header();
}
inline void Detection3DArray::unsafe_arena_set_allocated_header(
    ::apollo::common::Header* header) {
  if (GetArenaForAllocation() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(_impl_.header_);
  }
  _impl_.header_ = header;
  if (header) {
    
  } else {
    
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:apollo.contrib.lgsvl_msgs.Detection3DArray.header)
}
inline ::apollo::common::Header* Detection3DArray::release_header() {
  
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
inline ::apollo::common::Header* Detection3DArray::unsafe_arena_release_header() {
  // @@protoc_insertion_point(field_release:apollo.contrib.lgsvl_msgs.Detection3DArray.header)
  
  ::apollo::common::Header* temp = _impl_.header_;
  _impl_.header_ = nullptr;
  return temp;
}
inline ::apollo::common::Header* Detection3DArray::_internal_mutable_header() {
  
  if (_impl_.header_ == nullptr) {
    auto* p = CreateMaybeMessage<::apollo::common::Header>(GetArenaForAllocation());
    _impl_.header_ = p;
  }
  return _impl_.header_;
}
inline ::apollo::common::Header* Detection3DArray::mutable_header() {
  ::apollo::common::Header* _msg = _internal_mutable_header();
  // @@protoc_insertion_point(field_mutable:apollo.contrib.lgsvl_msgs.Detection3DArray.header)
  return _msg;
}
inline void Detection3DArray::set_allocated_header(::apollo::common::Header* header) {
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
    
  } else {
    
  }
  _impl_.header_ = header;
  // @@protoc_insertion_point(field_set_allocated:apollo.contrib.lgsvl_msgs.Detection3DArray.header)
}

// repeated .apollo.contrib.lgsvl_msgs.Detection3D detections = 2;
inline int Detection3DArray::_internal_detections_size() const {
  return _impl_.detections_.size();
}
inline int Detection3DArray::detections_size() const {
  return _internal_detections_size();
}
inline ::apollo::contrib::lgsvl_msgs::Detection3D* Detection3DArray::mutable_detections(int index) {
  // @@protoc_insertion_point(field_mutable:apollo.contrib.lgsvl_msgs.Detection3DArray.detections)
  return _impl_.detections_.Mutable(index);
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::apollo::contrib::lgsvl_msgs::Detection3D >*
Detection3DArray::mutable_detections() {
  // @@protoc_insertion_point(field_mutable_list:apollo.contrib.lgsvl_msgs.Detection3DArray.detections)
  return &_impl_.detections_;
}
inline const ::apollo::contrib::lgsvl_msgs::Detection3D& Detection3DArray::_internal_detections(int index) const {
  return _impl_.detections_.Get(index);
}
inline const ::apollo::contrib::lgsvl_msgs::Detection3D& Detection3DArray::detections(int index) const {
  // @@protoc_insertion_point(field_get:apollo.contrib.lgsvl_msgs.Detection3DArray.detections)
  return _internal_detections(index);
}
inline ::apollo::contrib::lgsvl_msgs::Detection3D* Detection3DArray::_internal_add_detections() {
  return _impl_.detections_.Add();
}
inline ::apollo::contrib::lgsvl_msgs::Detection3D* Detection3DArray::add_detections() {
  ::apollo::contrib::lgsvl_msgs::Detection3D* _add = _internal_add_detections();
  // @@protoc_insertion_point(field_add:apollo.contrib.lgsvl_msgs.Detection3DArray.detections)
  return _add;
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::apollo::contrib::lgsvl_msgs::Detection3D >&
Detection3DArray::detections() const {
  // @@protoc_insertion_point(field_list:apollo.contrib.lgsvl_msgs.Detection3DArray.detections)
  return _impl_.detections_;
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace lgsvl_msgs
}  // namespace contrib
}  // namespace apollo

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_modules_2fcontrib_2flgsvl_5fmsgs_2fproto_2fdetection3darray_2eproto
