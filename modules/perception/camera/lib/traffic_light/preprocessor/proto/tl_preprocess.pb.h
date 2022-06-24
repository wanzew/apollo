// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/perception/camera/lib/traffic_light/preprocessor/proto/tl_preprocess.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2fpreprocessor_2fproto_2ftl_5fpreprocess_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2fpreprocessor_2fproto_2ftl_5fpreprocess_2eproto

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
#include <google/protobuf/repeated_field.h>  // IWYU pragma: export
#include <google/protobuf/extension_set.h>  // IWYU pragma: export
#include <google/protobuf/generated_enum_reflection.h>
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2fpreprocessor_2fproto_2ftl_5fpreprocess_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2fpreprocessor_2fproto_2ftl_5fpreprocess_2eproto {
  static const uint32_t offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2fpreprocessor_2fproto_2ftl_5fpreprocess_2eproto;
PROTOBUF_NAMESPACE_OPEN
PROTOBUF_NAMESPACE_CLOSE
namespace apollo {
namespace perception {
namespace camera {
namespace traffic_light {
namespace preprocess {

enum CameraID : int {
  LONG_FOCUS = 0,
  NARROW_FOCUS = 1,
  SHORT_FOCUS = 2,
  WIDE_FOCUS = 3,
  UNKNOWN = 4,
  CAMERA_ID_COUNT = 5
};
bool CameraID_IsValid(int value);
constexpr CameraID CameraID_MIN = LONG_FOCUS;
constexpr CameraID CameraID_MAX = CAMERA_ID_COUNT;
constexpr int CameraID_ARRAYSIZE = CameraID_MAX + 1;

const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* CameraID_descriptor();
template<typename T>
inline const std::string& CameraID_Name(T enum_t_value) {
  static_assert(::std::is_same<T, CameraID>::value ||
    ::std::is_integral<T>::value,
    "Incorrect type passed to function CameraID_Name.");
  return ::PROTOBUF_NAMESPACE_ID::internal::NameOfEnum(
    CameraID_descriptor(), enum_t_value);
}
inline bool CameraID_Parse(
    ::PROTOBUF_NAMESPACE_ID::ConstStringParam name, CameraID* value) {
  return ::PROTOBUF_NAMESPACE_ID::internal::ParseNamedEnum<CameraID>(
    CameraID_descriptor(), name, value);
}
// ===================================================================


// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace preprocess
}  // namespace traffic_light
}  // namespace camera
}  // namespace perception
}  // namespace apollo

PROTOBUF_NAMESPACE_OPEN

template <> struct is_proto_enum< ::apollo::perception::camera::traffic_light::preprocess::CameraID> : ::std::true_type {};
template <>
inline const EnumDescriptor* GetEnumDescriptor< ::apollo::perception::camera::traffic_light::preprocess::CameraID>() {
  return ::apollo::perception::camera::traffic_light::preprocess::CameraID_descriptor();
}

PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_modules_2fperception_2fcamera_2flib_2ftraffic_5flight_2fpreprocessor_2fproto_2ftl_5fpreprocess_2eproto
