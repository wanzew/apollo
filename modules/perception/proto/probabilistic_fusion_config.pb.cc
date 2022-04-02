// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/perception/proto/probabilistic_fusion_config.proto

#include "modules/perception/proto/probabilistic_fusion_config.pb.h"

#include <algorithm>

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/extension_set.h>
#include <google/protobuf/wire_format_lite.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>

PROTOBUF_PRAGMA_INIT_SEG

namespace _pb = ::PROTOBUF_NAMESPACE_ID;
namespace _pbi = _pb::internal;

namespace apollo {
namespace perception {
namespace probabilistic_fusion_config {
PROTOBUF_CONSTEXPR ModelConfigs::ModelConfigs(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_._has_bits_)*/{}
  , /*decltype(_impl_._cached_size_)*/{}
  , /*decltype(_impl_.name_)*/{nullptr, ::_pbi::ConstantInitialized{}}
  , /*decltype(_impl_.version_)*/{nullptr, ::_pbi::ConstantInitialized{}}
  , /*decltype(_impl_.match_method_)*/{nullptr, ::_pbi::ConstantInitialized{}}
  , /*decltype(_impl_.publish_sensor_)*/{nullptr, ::_pbi::ConstantInitialized{}}
  , /*decltype(_impl_.max_camera_invisible_period_)*/0.25f
  , /*decltype(_impl_.max_match_distance_)*/4
  , /*decltype(_impl_.max_lidar_invisible_period_)*/0.25f
  , /*decltype(_impl_.max_radar_invisible_period_)*/0.25f
  , /*decltype(_impl_.max_radar_confident_angle_)*/30
  , /*decltype(_impl_.min_radar_confident_distance_)*/40
  , /*decltype(_impl_.publish_if_has_lidar_)*/true
  , /*decltype(_impl_.publish_if_has_radar_)*/true
  , /*decltype(_impl_.use_radar_)*/true
  , /*decltype(_impl_.use_lidar_)*/true} {}
struct ModelConfigsDefaultTypeInternal {
  PROTOBUF_CONSTEXPR ModelConfigsDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~ModelConfigsDefaultTypeInternal() {}
  union {
    ModelConfigs _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 ModelConfigsDefaultTypeInternal _ModelConfigs_default_instance_;
}  // namespace probabilistic_fusion_config
}  // namespace perception
}  // namespace apollo
static ::_pb::Metadata file_level_metadata_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto[1];
static constexpr ::_pb::EnumDescriptor const** file_level_enum_descriptors_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto = nullptr;
static constexpr ::_pb::ServiceDescriptor const** file_level_service_descriptors_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto = nullptr;

const uint32_t TableStruct_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_._has_bits_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.name_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.version_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.match_method_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.max_match_distance_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.max_lidar_invisible_period_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.max_radar_invisible_period_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.max_radar_confident_angle_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.min_radar_confident_distance_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.publish_if_has_lidar_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.publish_if_has_radar_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.publish_sensor_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.use_radar_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.use_lidar_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::probabilistic_fusion_config::ModelConfigs, _impl_.max_camera_invisible_period_),
  0,
  1,
  2,
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  3,
  12,
  13,
  4,
};
static const ::_pbi::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 20, -1, sizeof(::apollo::perception::probabilistic_fusion_config::ModelConfigs)},
};

static const ::_pb::Message* const file_default_instances[] = {
  &::apollo::perception::probabilistic_fusion_config::_ModelConfigs_default_instance_._instance,
};

const char descriptor_table_protodef_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n:modules/perception/proto/probabilistic"
  "_fusion_config.proto\022-apollo.perception."
  "probabilistic_fusion_config\"\371\003\n\014ModelCon"
  "figs\022!\n\004name\030\001 \001(\t:\023ProbabilisticFusion\022"
  "\026\n\007version\030\002 \001(\t:\0051.0.0\022 \n\014match_method\030"
  "\003 \001(\t:\nhm_matcher\022\035\n\022max_match_distance\030"
  "\004 \001(\002:\0014\022(\n\032max_lidar_invisible_period\030\005"
  " \001(\002:\0040.25\022(\n\032max_radar_invisible_period"
  "\030\006 \001(\002:\0040.25\022%\n\031max_radar_confident_angl"
  "e\030\007 \001(\002:\00230\022(\n\034min_radar_confident_dista"
  "nce\030\010 \001(\002:\00240\022\"\n\024publish_if_has_lidar\030\t "
  "\001(\010:\004true\022\"\n\024publish_if_has_radar\030\n \001(\010:"
  "\004true\022#\n\016publish_sensor\030\013 \001(\t:\013velodyne_"
  "64\022\027\n\tuse_radar\030\014 \001(\010:\004true\022\027\n\tuse_lidar"
  "\030\r \001(\010:\004true\022)\n\033max_camera_invisible_per"
  "iod\030\016 \001(\002:\0040.25"
  ;
static ::_pbi::once_flag descriptor_table_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto_once;
const ::_pbi::DescriptorTable descriptor_table_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto = {
    false, false, 615, descriptor_table_protodef_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto,
    "modules/perception/proto/probabilistic_fusion_config.proto",
    &descriptor_table_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto_once, nullptr, 0, 1,
    schemas, file_default_instances, TableStruct_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto::offsets,
    file_level_metadata_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto, file_level_enum_descriptors_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto,
    file_level_service_descriptors_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto,
};
PROTOBUF_ATTRIBUTE_WEAK const ::_pbi::DescriptorTable* descriptor_table_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto_getter() {
  return &descriptor_table_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto;
}

// Force running AddDescriptors() at dynamic initialization time.
PROTOBUF_ATTRIBUTE_INIT_PRIORITY2 static ::_pbi::AddDescriptorsRunner dynamic_init_dummy_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto(&descriptor_table_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto);
namespace apollo {
namespace perception {
namespace probabilistic_fusion_config {

// ===================================================================

class ModelConfigs::_Internal {
 public:
  using HasBits = decltype(std::declval<ModelConfigs>()._impl_._has_bits_);
  static void set_has_name(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_version(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static void set_has_match_method(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static void set_has_max_match_distance(HasBits* has_bits) {
    (*has_bits)[0] |= 32u;
  }
  static void set_has_max_lidar_invisible_period(HasBits* has_bits) {
    (*has_bits)[0] |= 64u;
  }
  static void set_has_max_radar_invisible_period(HasBits* has_bits) {
    (*has_bits)[0] |= 128u;
  }
  static void set_has_max_radar_confident_angle(HasBits* has_bits) {
    (*has_bits)[0] |= 256u;
  }
  static void set_has_min_radar_confident_distance(HasBits* has_bits) {
    (*has_bits)[0] |= 512u;
  }
  static void set_has_publish_if_has_lidar(HasBits* has_bits) {
    (*has_bits)[0] |= 1024u;
  }
  static void set_has_publish_if_has_radar(HasBits* has_bits) {
    (*has_bits)[0] |= 2048u;
  }
  static void set_has_publish_sensor(HasBits* has_bits) {
    (*has_bits)[0] |= 8u;
  }
  static void set_has_use_radar(HasBits* has_bits) {
    (*has_bits)[0] |= 4096u;
  }
  static void set_has_use_lidar(HasBits* has_bits) {
    (*has_bits)[0] |= 8192u;
  }
  static void set_has_max_camera_invisible_period(HasBits* has_bits) {
    (*has_bits)[0] |= 16u;
  }
};

const ::PROTOBUF_NAMESPACE_ID::internal::LazyString ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_name_{{{"ProbabilisticFusion", 19}}, {nullptr}};
const ::PROTOBUF_NAMESPACE_ID::internal::LazyString ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_version_{{{"1.0.0", 5}}, {nullptr}};
const ::PROTOBUF_NAMESPACE_ID::internal::LazyString ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_match_method_{{{"hm_matcher", 10}}, {nullptr}};
const ::PROTOBUF_NAMESPACE_ID::internal::LazyString ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_publish_sensor_{{{"velodyne_64", 11}}, {nullptr}};
ModelConfigs::ModelConfigs(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:apollo.perception.probabilistic_fusion_config.ModelConfigs)
}
ModelConfigs::ModelConfigs(const ModelConfigs& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){from._impl_._has_bits_}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.name_){}
    , decltype(_impl_.version_){}
    , decltype(_impl_.match_method_){}
    , decltype(_impl_.publish_sensor_){}
    , decltype(_impl_.max_camera_invisible_period_){}
    , decltype(_impl_.max_match_distance_){}
    , decltype(_impl_.max_lidar_invisible_period_){}
    , decltype(_impl_.max_radar_invisible_period_){}
    , decltype(_impl_.max_radar_confident_angle_){}
    , decltype(_impl_.min_radar_confident_distance_){}
    , decltype(_impl_.publish_if_has_lidar_){}
    , decltype(_impl_.publish_if_has_radar_){}
    , decltype(_impl_.use_radar_){}
    , decltype(_impl_.use_lidar_){}};

  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  _impl_.name_.InitDefault();
  if (from._internal_has_name()) {
    _impl_.name_.Set(from._internal_name(), 
      GetArenaForAllocation());
  }
  _impl_.version_.InitDefault();
  if (from._internal_has_version()) {
    _impl_.version_.Set(from._internal_version(), 
      GetArenaForAllocation());
  }
  _impl_.match_method_.InitDefault();
  if (from._internal_has_match_method()) {
    _impl_.match_method_.Set(from._internal_match_method(), 
      GetArenaForAllocation());
  }
  _impl_.publish_sensor_.InitDefault();
  if (from._internal_has_publish_sensor()) {
    _impl_.publish_sensor_.Set(from._internal_publish_sensor(), 
      GetArenaForAllocation());
  }
  ::memcpy(&_impl_.max_camera_invisible_period_, &from._impl_.max_camera_invisible_period_,
    static_cast<size_t>(reinterpret_cast<char*>(&_impl_.use_lidar_) -
    reinterpret_cast<char*>(&_impl_.max_camera_invisible_period_)) + sizeof(_impl_.use_lidar_));
  // @@protoc_insertion_point(copy_constructor:apollo.perception.probabilistic_fusion_config.ModelConfigs)
}

inline void ModelConfigs::SharedCtor(
    ::_pb::Arena* arena, bool is_message_owned) {
  (void)arena;
  (void)is_message_owned;
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.name_){}
    , decltype(_impl_.version_){}
    , decltype(_impl_.match_method_){}
    , decltype(_impl_.publish_sensor_){}
    , decltype(_impl_.max_camera_invisible_period_){0.25f}
    , decltype(_impl_.max_match_distance_){4}
    , decltype(_impl_.max_lidar_invisible_period_){0.25f}
    , decltype(_impl_.max_radar_invisible_period_){0.25f}
    , decltype(_impl_.max_radar_confident_angle_){30}
    , decltype(_impl_.min_radar_confident_distance_){40}
    , decltype(_impl_.publish_if_has_lidar_){true}
    , decltype(_impl_.publish_if_has_radar_){true}
    , decltype(_impl_.use_radar_){true}
    , decltype(_impl_.use_lidar_){true}
  };
  _impl_.name_.InitDefault();
  _impl_.version_.InitDefault();
  _impl_.match_method_.InitDefault();
  _impl_.publish_sensor_.InitDefault();
}

ModelConfigs::~ModelConfigs() {
  // @@protoc_insertion_point(destructor:apollo.perception.probabilistic_fusion_config.ModelConfigs)
  if (auto *arena = _internal_metadata_.DeleteReturnArena<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>()) {
  (void)arena;
    return;
  }
  SharedDtor();
}

inline void ModelConfigs::SharedDtor() {
  GOOGLE_DCHECK(GetArenaForAllocation() == nullptr);
  _impl_.name_.Destroy();
  _impl_.version_.Destroy();
  _impl_.match_method_.Destroy();
  _impl_.publish_sensor_.Destroy();
}

void ModelConfigs::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void ModelConfigs::Clear() {
// @@protoc_insertion_point(message_clear_start:apollo.perception.probabilistic_fusion_config.ModelConfigs)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x000000ffu) {
    if (cached_has_bits & 0x00000001u) {
      _impl_.name_.ClearToDefault(::apollo::perception::probabilistic_fusion_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_name_, GetArenaForAllocation());
       }
    if (cached_has_bits & 0x00000002u) {
      _impl_.version_.ClearToDefault(::apollo::perception::probabilistic_fusion_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_version_, GetArenaForAllocation());
       }
    if (cached_has_bits & 0x00000004u) {
      _impl_.match_method_.ClearToDefault(::apollo::perception::probabilistic_fusion_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_match_method_, GetArenaForAllocation());
       }
    if (cached_has_bits & 0x00000008u) {
      _impl_.publish_sensor_.ClearToDefault(::apollo::perception::probabilistic_fusion_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_publish_sensor_, GetArenaForAllocation());
       }
    _impl_.max_camera_invisible_period_ = 0.25f;
    _impl_.max_match_distance_ = 4;
    _impl_.max_lidar_invisible_period_ = 0.25f;
    _impl_.max_radar_invisible_period_ = 0.25f;
  }
  if (cached_has_bits & 0x00003f00u) {
    _impl_.max_radar_confident_angle_ = 30;
    _impl_.min_radar_confident_distance_ = 40;
    _impl_.publish_if_has_lidar_ = true;
    _impl_.publish_if_has_radar_ = true;
    _impl_.use_radar_ = true;
    _impl_.use_lidar_ = true;
  }
  _impl_._has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* ModelConfigs::_InternalParse(const char* ptr, ::_pbi::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  while (!ctx->Done(&ptr)) {
    uint32_t tag;
    ptr = ::_pbi::ReadTag(ptr, &tag);
    switch (tag >> 3) {
      // optional string name = 1 [default = "ProbabilisticFusion"];
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 10)) {
          auto str = _internal_mutable_name();
          ptr = ::_pbi::InlineGreedyStringParser(str, ptr, ctx);
          CHK_(ptr);
          #ifndef NDEBUG
          ::_pbi::VerifyUTF8(str, "apollo.perception.probabilistic_fusion_config.ModelConfigs.name");
          #endif  // !NDEBUG
        } else
          goto handle_unusual;
        continue;
      // optional string version = 2 [default = "1.0.0"];
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 18)) {
          auto str = _internal_mutable_version();
          ptr = ::_pbi::InlineGreedyStringParser(str, ptr, ctx);
          CHK_(ptr);
          #ifndef NDEBUG
          ::_pbi::VerifyUTF8(str, "apollo.perception.probabilistic_fusion_config.ModelConfigs.version");
          #endif  // !NDEBUG
        } else
          goto handle_unusual;
        continue;
      // optional string match_method = 3 [default = "hm_matcher"];
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 26)) {
          auto str = _internal_mutable_match_method();
          ptr = ::_pbi::InlineGreedyStringParser(str, ptr, ctx);
          CHK_(ptr);
          #ifndef NDEBUG
          ::_pbi::VerifyUTF8(str, "apollo.perception.probabilistic_fusion_config.ModelConfigs.match_method");
          #endif  // !NDEBUG
        } else
          goto handle_unusual;
        continue;
      // optional float max_match_distance = 4 [default = 4];
      case 4:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 37)) {
          _Internal::set_has_max_match_distance(&has_bits);
          _impl_.max_match_distance_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<float>(ptr);
          ptr += sizeof(float);
        } else
          goto handle_unusual;
        continue;
      // optional float max_lidar_invisible_period = 5 [default = 0.25];
      case 5:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 45)) {
          _Internal::set_has_max_lidar_invisible_period(&has_bits);
          _impl_.max_lidar_invisible_period_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<float>(ptr);
          ptr += sizeof(float);
        } else
          goto handle_unusual;
        continue;
      // optional float max_radar_invisible_period = 6 [default = 0.25];
      case 6:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 53)) {
          _Internal::set_has_max_radar_invisible_period(&has_bits);
          _impl_.max_radar_invisible_period_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<float>(ptr);
          ptr += sizeof(float);
        } else
          goto handle_unusual;
        continue;
      // optional float max_radar_confident_angle = 7 [default = 30];
      case 7:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 61)) {
          _Internal::set_has_max_radar_confident_angle(&has_bits);
          _impl_.max_radar_confident_angle_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<float>(ptr);
          ptr += sizeof(float);
        } else
          goto handle_unusual;
        continue;
      // optional float min_radar_confident_distance = 8 [default = 40];
      case 8:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 69)) {
          _Internal::set_has_min_radar_confident_distance(&has_bits);
          _impl_.min_radar_confident_distance_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<float>(ptr);
          ptr += sizeof(float);
        } else
          goto handle_unusual;
        continue;
      // optional bool publish_if_has_lidar = 9 [default = true];
      case 9:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 72)) {
          _Internal::set_has_publish_if_has_lidar(&has_bits);
          _impl_.publish_if_has_lidar_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional bool publish_if_has_radar = 10 [default = true];
      case 10:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 80)) {
          _Internal::set_has_publish_if_has_radar(&has_bits);
          _impl_.publish_if_has_radar_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional string publish_sensor = 11 [default = "velodyne_64"];
      case 11:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 90)) {
          auto str = _internal_mutable_publish_sensor();
          ptr = ::_pbi::InlineGreedyStringParser(str, ptr, ctx);
          CHK_(ptr);
          #ifndef NDEBUG
          ::_pbi::VerifyUTF8(str, "apollo.perception.probabilistic_fusion_config.ModelConfigs.publish_sensor");
          #endif  // !NDEBUG
        } else
          goto handle_unusual;
        continue;
      // optional bool use_radar = 12 [default = true];
      case 12:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 96)) {
          _Internal::set_has_use_radar(&has_bits);
          _impl_.use_radar_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional bool use_lidar = 13 [default = true];
      case 13:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 104)) {
          _Internal::set_has_use_lidar(&has_bits);
          _impl_.use_lidar_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // optional float max_camera_invisible_period = 14 [default = 0.25];
      case 14:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 117)) {
          _Internal::set_has_max_camera_invisible_period(&has_bits);
          _impl_.max_camera_invisible_period_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<float>(ptr);
          ptr += sizeof(float);
        } else
          goto handle_unusual;
        continue;
      default:
        goto handle_unusual;
    }  // switch
  handle_unusual:
    if ((tag == 0) || ((tag & 7) == 4)) {
      CHK_(ptr);
      ctx->SetLastTag(tag);
      goto message_done;
    }
    ptr = UnknownFieldParse(
        tag,
        _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(),
        ptr, ctx);
    CHK_(ptr != nullptr);
  }  // while
message_done:
  _impl_._has_bits_.Or(has_bits);
  return ptr;
failure:
  ptr = nullptr;
  goto message_done;
#undef CHK_
}

uint8_t* ModelConfigs::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:apollo.perception.probabilistic_fusion_config.ModelConfigs)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  // optional string name = 1 [default = "ProbabilisticFusion"];
  if (cached_has_bits & 0x00000001u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_name().data(), static_cast<int>(this->_internal_name().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "apollo.perception.probabilistic_fusion_config.ModelConfigs.name");
    target = stream->WriteStringMaybeAliased(
        1, this->_internal_name(), target);
  }

  // optional string version = 2 [default = "1.0.0"];
  if (cached_has_bits & 0x00000002u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_version().data(), static_cast<int>(this->_internal_version().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "apollo.perception.probabilistic_fusion_config.ModelConfigs.version");
    target = stream->WriteStringMaybeAliased(
        2, this->_internal_version(), target);
  }

  // optional string match_method = 3 [default = "hm_matcher"];
  if (cached_has_bits & 0x00000004u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_match_method().data(), static_cast<int>(this->_internal_match_method().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "apollo.perception.probabilistic_fusion_config.ModelConfigs.match_method");
    target = stream->WriteStringMaybeAliased(
        3, this->_internal_match_method(), target);
  }

  // optional float max_match_distance = 4 [default = 4];
  if (cached_has_bits & 0x00000020u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteFloatToArray(4, this->_internal_max_match_distance(), target);
  }

  // optional float max_lidar_invisible_period = 5 [default = 0.25];
  if (cached_has_bits & 0x00000040u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteFloatToArray(5, this->_internal_max_lidar_invisible_period(), target);
  }

  // optional float max_radar_invisible_period = 6 [default = 0.25];
  if (cached_has_bits & 0x00000080u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteFloatToArray(6, this->_internal_max_radar_invisible_period(), target);
  }

  // optional float max_radar_confident_angle = 7 [default = 30];
  if (cached_has_bits & 0x00000100u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteFloatToArray(7, this->_internal_max_radar_confident_angle(), target);
  }

  // optional float min_radar_confident_distance = 8 [default = 40];
  if (cached_has_bits & 0x00000200u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteFloatToArray(8, this->_internal_min_radar_confident_distance(), target);
  }

  // optional bool publish_if_has_lidar = 9 [default = true];
  if (cached_has_bits & 0x00000400u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteBoolToArray(9, this->_internal_publish_if_has_lidar(), target);
  }

  // optional bool publish_if_has_radar = 10 [default = true];
  if (cached_has_bits & 0x00000800u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteBoolToArray(10, this->_internal_publish_if_has_radar(), target);
  }

  // optional string publish_sensor = 11 [default = "velodyne_64"];
  if (cached_has_bits & 0x00000008u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_publish_sensor().data(), static_cast<int>(this->_internal_publish_sensor().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "apollo.perception.probabilistic_fusion_config.ModelConfigs.publish_sensor");
    target = stream->WriteStringMaybeAliased(
        11, this->_internal_publish_sensor(), target);
  }

  // optional bool use_radar = 12 [default = true];
  if (cached_has_bits & 0x00001000u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteBoolToArray(12, this->_internal_use_radar(), target);
  }

  // optional bool use_lidar = 13 [default = true];
  if (cached_has_bits & 0x00002000u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteBoolToArray(13, this->_internal_use_lidar(), target);
  }

  // optional float max_camera_invisible_period = 14 [default = 0.25];
  if (cached_has_bits & 0x00000010u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteFloatToArray(14, this->_internal_max_camera_invisible_period(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:apollo.perception.probabilistic_fusion_config.ModelConfigs)
  return target;
}

size_t ModelConfigs::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:apollo.perception.probabilistic_fusion_config.ModelConfigs)
  size_t total_size = 0;

  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x000000ffu) {
    // optional string name = 1 [default = "ProbabilisticFusion"];
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
          this->_internal_name());
    }

    // optional string version = 2 [default = "1.0.0"];
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
          this->_internal_version());
    }

    // optional string match_method = 3 [default = "hm_matcher"];
    if (cached_has_bits & 0x00000004u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
          this->_internal_match_method());
    }

    // optional string publish_sensor = 11 [default = "velodyne_64"];
    if (cached_has_bits & 0x00000008u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
          this->_internal_publish_sensor());
    }

    // optional float max_camera_invisible_period = 14 [default = 0.25];
    if (cached_has_bits & 0x00000010u) {
      total_size += 1 + 4;
    }

    // optional float max_match_distance = 4 [default = 4];
    if (cached_has_bits & 0x00000020u) {
      total_size += 1 + 4;
    }

    // optional float max_lidar_invisible_period = 5 [default = 0.25];
    if (cached_has_bits & 0x00000040u) {
      total_size += 1 + 4;
    }

    // optional float max_radar_invisible_period = 6 [default = 0.25];
    if (cached_has_bits & 0x00000080u) {
      total_size += 1 + 4;
    }

  }
  if (cached_has_bits & 0x00003f00u) {
    // optional float max_radar_confident_angle = 7 [default = 30];
    if (cached_has_bits & 0x00000100u) {
      total_size += 1 + 4;
    }

    // optional float min_radar_confident_distance = 8 [default = 40];
    if (cached_has_bits & 0x00000200u) {
      total_size += 1 + 4;
    }

    // optional bool publish_if_has_lidar = 9 [default = true];
    if (cached_has_bits & 0x00000400u) {
      total_size += 1 + 1;
    }

    // optional bool publish_if_has_radar = 10 [default = true];
    if (cached_has_bits & 0x00000800u) {
      total_size += 1 + 1;
    }

    // optional bool use_radar = 12 [default = true];
    if (cached_has_bits & 0x00001000u) {
      total_size += 1 + 1;
    }

    // optional bool use_lidar = 13 [default = true];
    if (cached_has_bits & 0x00002000u) {
      total_size += 1 + 1;
    }

  }
  return MaybeComputeUnknownFieldsSize(total_size, &_impl_._cached_size_);
}

const ::PROTOBUF_NAMESPACE_ID::Message::ClassData ModelConfigs::_class_data_ = {
    ::PROTOBUF_NAMESPACE_ID::Message::CopyWithSizeCheck,
    ModelConfigs::MergeImpl
};
const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*ModelConfigs::GetClassData() const { return &_class_data_; }

void ModelConfigs::MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to,
                      const ::PROTOBUF_NAMESPACE_ID::Message& from) {
  static_cast<ModelConfigs *>(to)->MergeFrom(
      static_cast<const ModelConfigs &>(from));
}


void ModelConfigs::MergeFrom(const ModelConfigs& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:apollo.perception.probabilistic_fusion_config.ModelConfigs)
  GOOGLE_DCHECK_NE(&from, this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._impl_._has_bits_[0];
  if (cached_has_bits & 0x000000ffu) {
    if (cached_has_bits & 0x00000001u) {
      _internal_set_name(from._internal_name());
    }
    if (cached_has_bits & 0x00000002u) {
      _internal_set_version(from._internal_version());
    }
    if (cached_has_bits & 0x00000004u) {
      _internal_set_match_method(from._internal_match_method());
    }
    if (cached_has_bits & 0x00000008u) {
      _internal_set_publish_sensor(from._internal_publish_sensor());
    }
    if (cached_has_bits & 0x00000010u) {
      _impl_.max_camera_invisible_period_ = from._impl_.max_camera_invisible_period_;
    }
    if (cached_has_bits & 0x00000020u) {
      _impl_.max_match_distance_ = from._impl_.max_match_distance_;
    }
    if (cached_has_bits & 0x00000040u) {
      _impl_.max_lidar_invisible_period_ = from._impl_.max_lidar_invisible_period_;
    }
    if (cached_has_bits & 0x00000080u) {
      _impl_.max_radar_invisible_period_ = from._impl_.max_radar_invisible_period_;
    }
    _impl_._has_bits_[0] |= cached_has_bits;
  }
  if (cached_has_bits & 0x00003f00u) {
    if (cached_has_bits & 0x00000100u) {
      _impl_.max_radar_confident_angle_ = from._impl_.max_radar_confident_angle_;
    }
    if (cached_has_bits & 0x00000200u) {
      _impl_.min_radar_confident_distance_ = from._impl_.min_radar_confident_distance_;
    }
    if (cached_has_bits & 0x00000400u) {
      _impl_.publish_if_has_lidar_ = from._impl_.publish_if_has_lidar_;
    }
    if (cached_has_bits & 0x00000800u) {
      _impl_.publish_if_has_radar_ = from._impl_.publish_if_has_radar_;
    }
    if (cached_has_bits & 0x00001000u) {
      _impl_.use_radar_ = from._impl_.use_radar_;
    }
    if (cached_has_bits & 0x00002000u) {
      _impl_.use_lidar_ = from._impl_.use_lidar_;
    }
    _impl_._has_bits_[0] |= cached_has_bits;
  }
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void ModelConfigs::CopyFrom(const ModelConfigs& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:apollo.perception.probabilistic_fusion_config.ModelConfigs)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool ModelConfigs::IsInitialized() const {
  return true;
}

void ModelConfigs::InternalSwap(ModelConfigs* other) {
  using std::swap;
  auto* lhs_arena = GetArenaForAllocation();
  auto* rhs_arena = other->GetArenaForAllocation();
  _internal_metadata_.InternalSwap(&other->_internal_metadata_);
  swap(_impl_._has_bits_[0], other->_impl_._has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr::InternalSwap(
      &_impl_.name_, lhs_arena,
      &other->_impl_.name_, rhs_arena
  );
  ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr::InternalSwap(
      &_impl_.version_, lhs_arena,
      &other->_impl_.version_, rhs_arena
  );
  ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr::InternalSwap(
      &_impl_.match_method_, lhs_arena,
      &other->_impl_.match_method_, rhs_arena
  );
  ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr::InternalSwap(
      &_impl_.publish_sensor_, lhs_arena,
      &other->_impl_.publish_sensor_, rhs_arena
  );
  swap(_impl_.max_camera_invisible_period_, other->_impl_.max_camera_invisible_period_);
  swap(_impl_.max_match_distance_, other->_impl_.max_match_distance_);
  swap(_impl_.max_lidar_invisible_period_, other->_impl_.max_lidar_invisible_period_);
  swap(_impl_.max_radar_invisible_period_, other->_impl_.max_radar_invisible_period_);
  swap(_impl_.max_radar_confident_angle_, other->_impl_.max_radar_confident_angle_);
  swap(_impl_.min_radar_confident_distance_, other->_impl_.min_radar_confident_distance_);
  swap(_impl_.publish_if_has_lidar_, other->_impl_.publish_if_has_lidar_);
  swap(_impl_.publish_if_has_radar_, other->_impl_.publish_if_has_radar_);
  swap(_impl_.use_radar_, other->_impl_.use_radar_);
  swap(_impl_.use_lidar_, other->_impl_.use_lidar_);
}

::PROTOBUF_NAMESPACE_ID::Metadata ModelConfigs::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto_getter, &descriptor_table_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto_once,
      file_level_metadata_modules_2fperception_2fproto_2fprobabilistic_5ffusion_5fconfig_2eproto[0]);
}

// @@protoc_insertion_point(namespace_scope)
}  // namespace probabilistic_fusion_config
}  // namespace perception
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::apollo::perception::probabilistic_fusion_config::ModelConfigs*
Arena::CreateMaybeMessage< ::apollo::perception::probabilistic_fusion_config::ModelConfigs >(Arena* arena) {
  return Arena::CreateMessageInternal< ::apollo::perception::probabilistic_fusion_config::ModelConfigs >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
