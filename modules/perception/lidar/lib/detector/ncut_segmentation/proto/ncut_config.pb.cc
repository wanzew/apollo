// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/perception/lidar/lib/detector/ncut_segmentation/proto/ncut_config.proto

#include "modules/perception/lidar/lib/detector/ncut_segmentation/proto/ncut_config.pb.h"

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
namespace lidar {
PROTOBUF_CONSTEXPR NCutConfig::NCutConfig(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_._has_bits_)*/{}
  , /*decltype(_impl_._cached_size_)*/{}
  , /*decltype(_impl_.param_file_)*/{nullptr, ::_pbi::ConstantInitialized{}}} {}
struct NCutConfigDefaultTypeInternal {
  PROTOBUF_CONSTEXPR NCutConfigDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~NCutConfigDefaultTypeInternal() {}
  union {
    NCutConfig _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 NCutConfigDefaultTypeInternal _NCutConfig_default_instance_;
}  // namespace lidar
}  // namespace perception
}  // namespace apollo
static ::_pb::Metadata file_level_metadata_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto[1];
static constexpr ::_pb::EnumDescriptor const** file_level_enum_descriptors_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto = nullptr;
static constexpr ::_pb::ServiceDescriptor const** file_level_service_descriptors_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto = nullptr;

const uint32_t TableStruct_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::apollo::perception::lidar::NCutConfig, _impl_._has_bits_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::lidar::NCutConfig, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::apollo::perception::lidar::NCutConfig, _impl_.param_file_),
  0,
};
static const ::_pbi::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 7, -1, sizeof(::apollo::perception::lidar::NCutConfig)},
};

static const ::_pb::Message* const file_default_instances[] = {
  &::apollo::perception::lidar::_NCutConfig_default_instance_._instance,
};

const char descriptor_table_protodef_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\nOmodules/perception/lidar/lib/detector/"
  "ncut_segmentation/proto/ncut_config.prot"
  "o\022\027apollo.perception.lidar\"\?\n\nNCutConfig"
  "\0221\n\nparam_file\030\001 \001(\t:\035./data/models/ncut"
  "/param.conf"
  ;
static ::_pbi::once_flag descriptor_table_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto_once;
const ::_pbi::DescriptorTable descriptor_table_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto = {
    false, false, 171, descriptor_table_protodef_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto,
    "modules/perception/lidar/lib/detector/ncut_segmentation/proto/ncut_config.proto",
    &descriptor_table_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto_once, nullptr, 0, 1,
    schemas, file_default_instances, TableStruct_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto::offsets,
    file_level_metadata_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto, file_level_enum_descriptors_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto,
    file_level_service_descriptors_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto,
};
PROTOBUF_ATTRIBUTE_WEAK const ::_pbi::DescriptorTable* descriptor_table_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto_getter() {
  return &descriptor_table_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto;
}

// Force running AddDescriptors() at dynamic initialization time.
PROTOBUF_ATTRIBUTE_INIT_PRIORITY2 static ::_pbi::AddDescriptorsRunner dynamic_init_dummy_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto(&descriptor_table_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto);
namespace apollo {
namespace perception {
namespace lidar {

// ===================================================================

class NCutConfig::_Internal {
 public:
  using HasBits = decltype(std::declval<NCutConfig>()._impl_._has_bits_);
  static void set_has_param_file(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
};

const ::PROTOBUF_NAMESPACE_ID::internal::LazyString NCutConfig::Impl_::_i_give_permission_to_break_this_code_default_param_file_{{{"./data/models/ncut/param.conf", 29}}, {nullptr}};
NCutConfig::NCutConfig(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:apollo.perception.lidar.NCutConfig)
}
NCutConfig::NCutConfig(const NCutConfig& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){from._impl_._has_bits_}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.param_file_){}};

  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  _impl_.param_file_.InitDefault();
  if (from._internal_has_param_file()) {
    _impl_.param_file_.Set(from._internal_param_file(), 
      GetArenaForAllocation());
  }
  // @@protoc_insertion_point(copy_constructor:apollo.perception.lidar.NCutConfig)
}

inline void NCutConfig::SharedCtor(
    ::_pb::Arena* arena, bool is_message_owned) {
  (void)arena;
  (void)is_message_owned;
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.param_file_){}
  };
  _impl_.param_file_.InitDefault();
}

NCutConfig::~NCutConfig() {
  // @@protoc_insertion_point(destructor:apollo.perception.lidar.NCutConfig)
  if (auto *arena = _internal_metadata_.DeleteReturnArena<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>()) {
  (void)arena;
    return;
  }
  SharedDtor();
}

inline void NCutConfig::SharedDtor() {
  GOOGLE_DCHECK(GetArenaForAllocation() == nullptr);
  _impl_.param_file_.Destroy();
}

void NCutConfig::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void NCutConfig::Clear() {
// @@protoc_insertion_point(message_clear_start:apollo.perception.lidar.NCutConfig)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x00000001u) {
    _impl_.param_file_.ClearToDefault(::apollo::perception::lidar::NCutConfig::Impl_::_i_give_permission_to_break_this_code_default_param_file_, GetArenaForAllocation());
     }
  _impl_._has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* NCutConfig::_InternalParse(const char* ptr, ::_pbi::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  while (!ctx->Done(&ptr)) {
    uint32_t tag;
    ptr = ::_pbi::ReadTag(ptr, &tag);
    switch (tag >> 3) {
      // optional string param_file = 1 [default = "./data/models/ncut/param.conf"];
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 10)) {
          auto str = _internal_mutable_param_file();
          ptr = ::_pbi::InlineGreedyStringParser(str, ptr, ctx);
          CHK_(ptr);
          #ifndef NDEBUG
          ::_pbi::VerifyUTF8(str, "apollo.perception.lidar.NCutConfig.param_file");
          #endif  // !NDEBUG
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

uint8_t* NCutConfig::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:apollo.perception.lidar.NCutConfig)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  // optional string param_file = 1 [default = "./data/models/ncut/param.conf"];
  if (cached_has_bits & 0x00000001u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_param_file().data(), static_cast<int>(this->_internal_param_file().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "apollo.perception.lidar.NCutConfig.param_file");
    target = stream->WriteStringMaybeAliased(
        1, this->_internal_param_file(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:apollo.perception.lidar.NCutConfig)
  return target;
}

size_t NCutConfig::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:apollo.perception.lidar.NCutConfig)
  size_t total_size = 0;

  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // optional string param_file = 1 [default = "./data/models/ncut/param.conf"];
  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x00000001u) {
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
        this->_internal_param_file());
  }

  return MaybeComputeUnknownFieldsSize(total_size, &_impl_._cached_size_);
}

const ::PROTOBUF_NAMESPACE_ID::Message::ClassData NCutConfig::_class_data_ = {
    ::PROTOBUF_NAMESPACE_ID::Message::CopyWithSizeCheck,
    NCutConfig::MergeImpl
};
const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*NCutConfig::GetClassData() const { return &_class_data_; }

void NCutConfig::MergeImpl(::PROTOBUF_NAMESPACE_ID::Message* to,
                      const ::PROTOBUF_NAMESPACE_ID::Message& from) {
  static_cast<NCutConfig *>(to)->MergeFrom(
      static_cast<const NCutConfig &>(from));
}


void NCutConfig::MergeFrom(const NCutConfig& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:apollo.perception.lidar.NCutConfig)
  GOOGLE_DCHECK_NE(&from, this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  if (from._internal_has_param_file()) {
    _internal_set_param_file(from._internal_param_file());
  }
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void NCutConfig::CopyFrom(const NCutConfig& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:apollo.perception.lidar.NCutConfig)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool NCutConfig::IsInitialized() const {
  return true;
}

void NCutConfig::InternalSwap(NCutConfig* other) {
  using std::swap;
  auto* lhs_arena = GetArenaForAllocation();
  auto* rhs_arena = other->GetArenaForAllocation();
  _internal_metadata_.InternalSwap(&other->_internal_metadata_);
  swap(_impl_._has_bits_[0], other->_impl_._has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr::InternalSwap(
      &_impl_.param_file_, lhs_arena,
      &other->_impl_.param_file_, rhs_arena
  );
}

::PROTOBUF_NAMESPACE_ID::Metadata NCutConfig::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto_getter, &descriptor_table_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto_once,
      file_level_metadata_modules_2fperception_2flidar_2flib_2fdetector_2fncut_5fsegmentation_2fproto_2fncut_5fconfig_2eproto[0]);
}

// @@protoc_insertion_point(namespace_scope)
}  // namespace lidar
}  // namespace perception
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::apollo::perception::lidar::NCutConfig*
Arena::CreateMaybeMessage< ::apollo::perception::lidar::NCutConfig >(Arena* arena) {
  return Arena::CreateMessageInternal< ::apollo::perception::lidar::NCutConfig >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
