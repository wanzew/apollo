// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: modules/perception/proto/sequence_type_fuser_config.proto

#include "modules/perception/proto/sequence_type_fuser_config.pb.h"

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
namespace sequence_type_fuser_config {
PROTOBUF_CONSTEXPR ModelConfigs::ModelConfigs(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_._has_bits_)*/{}
  , /*decltype(_impl_._cached_size_)*/{}
  , /*decltype(_impl_.name_)*/{nullptr, ::_pbi::ConstantInitialized{}}
  , /*decltype(_impl_.version_)*/{nullptr, ::_pbi::ConstantInitialized{}}
  , /*decltype(_impl_.classifiers_property_file_path_)*/{nullptr, ::_pbi::ConstantInitialized{}}
  , /*decltype(_impl_.transition_property_file_path_)*/{nullptr, ::_pbi::ConstantInitialized{}}
  , /*decltype(_impl_.temporal_window_)*/20} {}
struct ModelConfigsDefaultTypeInternal {
  PROTOBUF_CONSTEXPR ModelConfigsDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~ModelConfigsDefaultTypeInternal() {}
  union {
    ModelConfigs _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 ModelConfigsDefaultTypeInternal _ModelConfigs_default_instance_;
}  // namespace sequence_type_fuser_config
}  // namespace perception
}  // namespace apollo
static ::_pb::Metadata file_level_metadata_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto[1];
static constexpr ::_pb::EnumDescriptor const** file_level_enum_descriptors_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto = nullptr;
static constexpr ::_pb::ServiceDescriptor const** file_level_service_descriptors_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto = nullptr;

const uint32_t TableStruct_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::apollo::perception::sequence_type_fuser_config::ModelConfigs, _impl_._has_bits_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::sequence_type_fuser_config::ModelConfigs, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::apollo::perception::sequence_type_fuser_config::ModelConfigs, _impl_.name_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::sequence_type_fuser_config::ModelConfigs, _impl_.version_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::sequence_type_fuser_config::ModelConfigs, _impl_.temporal_window_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::sequence_type_fuser_config::ModelConfigs, _impl_.classifiers_property_file_path_),
  PROTOBUF_FIELD_OFFSET(::apollo::perception::sequence_type_fuser_config::ModelConfigs, _impl_.transition_property_file_path_),
  0,
  1,
  4,
  2,
  3,
};
static const ::_pbi::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 11, -1, sizeof(::apollo::perception::sequence_type_fuser_config::ModelConfigs)},
};

static const ::_pb::Message* const file_default_instances[] = {
  &::apollo::perception::sequence_type_fuser_config::_ModelConfigs_default_instance_._instance,
};

const char descriptor_table_protodef_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n9modules/perception/proto/sequence_type"
  "_fuser_config.proto\022,apollo.perception.s"
  "equence_type_fuser_config\"\270\002\n\014ModelConfi"
  "gs\022\037\n\004name\030\001 \001(\t:\021SequenceTypeFuser\022\026\n\007v"
  "ersion\030\002 \001(\t:\0051.1.0\022\033\n\017temporal_window\030\003"
  " \001(\002:\00220\022i\n\036classifiers_property_file_pa"
  "th\030\004 \001(\t:Amodules/perception/model/seque"
  "nce_type_fuser/classifiers.property\022g\n\035t"
  "ransition_property_file_path\030\005 \001(\t:@modu"
  "les/perception/model/sequence_type_fuser"
  "/transition.property"
  ;
static ::_pbi::once_flag descriptor_table_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto_once;
const ::_pbi::DescriptorTable descriptor_table_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto = {
    false, false, 420, descriptor_table_protodef_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto,
    "modules/perception/proto/sequence_type_fuser_config.proto",
    &descriptor_table_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto_once, nullptr, 0, 1,
    schemas, file_default_instances, TableStruct_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto::offsets,
    file_level_metadata_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto, file_level_enum_descriptors_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto,
    file_level_service_descriptors_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto,
};
PROTOBUF_ATTRIBUTE_WEAK const ::_pbi::DescriptorTable* descriptor_table_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto_getter() {
  return &descriptor_table_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto;
}

// Force running AddDescriptors() at dynamic initialization time.
PROTOBUF_ATTRIBUTE_INIT_PRIORITY2 static ::_pbi::AddDescriptorsRunner dynamic_init_dummy_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto(&descriptor_table_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto);
namespace apollo {
namespace perception {
namespace sequence_type_fuser_config {

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
  static void set_has_temporal_window(HasBits* has_bits) {
    (*has_bits)[0] |= 16u;
  }
  static void set_has_classifiers_property_file_path(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static void set_has_transition_property_file_path(HasBits* has_bits) {
    (*has_bits)[0] |= 8u;
  }
};

const ::PROTOBUF_NAMESPACE_ID::internal::LazyString ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_name_{{{"SequenceTypeFuser", 17}}, {nullptr}};
const ::PROTOBUF_NAMESPACE_ID::internal::LazyString ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_version_{{{"1.1.0", 5}}, {nullptr}};
const ::PROTOBUF_NAMESPACE_ID::internal::LazyString ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_classifiers_property_file_path_{{{"modules/perception/model/sequence_type_fuser/classifiers.property", 65}}, {nullptr}};
const ::PROTOBUF_NAMESPACE_ID::internal::LazyString ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_transition_property_file_path_{{{"modules/perception/model/sequence_type_fuser/transition.property", 64}}, {nullptr}};
ModelConfigs::ModelConfigs(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:apollo.perception.sequence_type_fuser_config.ModelConfigs)
}
ModelConfigs::ModelConfigs(const ModelConfigs& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  new (&_impl_) Impl_{
      decltype(_impl_._has_bits_){from._impl_._has_bits_}
    , /*decltype(_impl_._cached_size_)*/{}
    , decltype(_impl_.name_){}
    , decltype(_impl_.version_){}
    , decltype(_impl_.classifiers_property_file_path_){}
    , decltype(_impl_.transition_property_file_path_){}
    , decltype(_impl_.temporal_window_){}};

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
  _impl_.classifiers_property_file_path_.InitDefault();
  if (from._internal_has_classifiers_property_file_path()) {
    _impl_.classifiers_property_file_path_.Set(from._internal_classifiers_property_file_path(), 
      GetArenaForAllocation());
  }
  _impl_.transition_property_file_path_.InitDefault();
  if (from._internal_has_transition_property_file_path()) {
    _impl_.transition_property_file_path_.Set(from._internal_transition_property_file_path(), 
      GetArenaForAllocation());
  }
  _impl_.temporal_window_ = from._impl_.temporal_window_;
  // @@protoc_insertion_point(copy_constructor:apollo.perception.sequence_type_fuser_config.ModelConfigs)
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
    , decltype(_impl_.classifiers_property_file_path_){}
    , decltype(_impl_.transition_property_file_path_){}
    , decltype(_impl_.temporal_window_){20}
  };
  _impl_.name_.InitDefault();
  _impl_.version_.InitDefault();
  _impl_.classifiers_property_file_path_.InitDefault();
  _impl_.transition_property_file_path_.InitDefault();
}

ModelConfigs::~ModelConfigs() {
  // @@protoc_insertion_point(destructor:apollo.perception.sequence_type_fuser_config.ModelConfigs)
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
  _impl_.classifiers_property_file_path_.Destroy();
  _impl_.transition_property_file_path_.Destroy();
}

void ModelConfigs::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void ModelConfigs::Clear() {
// @@protoc_insertion_point(message_clear_start:apollo.perception.sequence_type_fuser_config.ModelConfigs)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    if (cached_has_bits & 0x00000001u) {
      _impl_.name_.ClearToDefault(::apollo::perception::sequence_type_fuser_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_name_, GetArenaForAllocation());
       }
    if (cached_has_bits & 0x00000002u) {
      _impl_.version_.ClearToDefault(::apollo::perception::sequence_type_fuser_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_version_, GetArenaForAllocation());
       }
    if (cached_has_bits & 0x00000004u) {
      _impl_.classifiers_property_file_path_.ClearToDefault(::apollo::perception::sequence_type_fuser_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_classifiers_property_file_path_, GetArenaForAllocation());
       }
    if (cached_has_bits & 0x00000008u) {
      _impl_.transition_property_file_path_.ClearToDefault(::apollo::perception::sequence_type_fuser_config::ModelConfigs::Impl_::_i_give_permission_to_break_this_code_default_transition_property_file_path_, GetArenaForAllocation());
       }
    _impl_.temporal_window_ = 20;
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
      // optional string name = 1 [default = "SequenceTypeFuser"];
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 10)) {
          auto str = _internal_mutable_name();
          ptr = ::_pbi::InlineGreedyStringParser(str, ptr, ctx);
          CHK_(ptr);
          #ifndef NDEBUG
          ::_pbi::VerifyUTF8(str, "apollo.perception.sequence_type_fuser_config.ModelConfigs.name");
          #endif  // !NDEBUG
        } else
          goto handle_unusual;
        continue;
      // optional string version = 2 [default = "1.1.0"];
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 18)) {
          auto str = _internal_mutable_version();
          ptr = ::_pbi::InlineGreedyStringParser(str, ptr, ctx);
          CHK_(ptr);
          #ifndef NDEBUG
          ::_pbi::VerifyUTF8(str, "apollo.perception.sequence_type_fuser_config.ModelConfigs.version");
          #endif  // !NDEBUG
        } else
          goto handle_unusual;
        continue;
      // optional float temporal_window = 3 [default = 20];
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 29)) {
          _Internal::set_has_temporal_window(&has_bits);
          _impl_.temporal_window_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<float>(ptr);
          ptr += sizeof(float);
        } else
          goto handle_unusual;
        continue;
      // optional string classifiers_property_file_path = 4 [default = "modules/perception/model/sequence_type_fuser/classifiers.property"];
      case 4:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 34)) {
          auto str = _internal_mutable_classifiers_property_file_path();
          ptr = ::_pbi::InlineGreedyStringParser(str, ptr, ctx);
          CHK_(ptr);
          #ifndef NDEBUG
          ::_pbi::VerifyUTF8(str, "apollo.perception.sequence_type_fuser_config.ModelConfigs.classifiers_property_file_path");
          #endif  // !NDEBUG
        } else
          goto handle_unusual;
        continue;
      // optional string transition_property_file_path = 5 [default = "modules/perception/model/sequence_type_fuser/transition.property"];
      case 5:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 42)) {
          auto str = _internal_mutable_transition_property_file_path();
          ptr = ::_pbi::InlineGreedyStringParser(str, ptr, ctx);
          CHK_(ptr);
          #ifndef NDEBUG
          ::_pbi::VerifyUTF8(str, "apollo.perception.sequence_type_fuser_config.ModelConfigs.transition_property_file_path");
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

uint8_t* ModelConfigs::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:apollo.perception.sequence_type_fuser_config.ModelConfigs)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  // optional string name = 1 [default = "SequenceTypeFuser"];
  if (cached_has_bits & 0x00000001u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_name().data(), static_cast<int>(this->_internal_name().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "apollo.perception.sequence_type_fuser_config.ModelConfigs.name");
    target = stream->WriteStringMaybeAliased(
        1, this->_internal_name(), target);
  }

  // optional string version = 2 [default = "1.1.0"];
  if (cached_has_bits & 0x00000002u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_version().data(), static_cast<int>(this->_internal_version().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "apollo.perception.sequence_type_fuser_config.ModelConfigs.version");
    target = stream->WriteStringMaybeAliased(
        2, this->_internal_version(), target);
  }

  // optional float temporal_window = 3 [default = 20];
  if (cached_has_bits & 0x00000010u) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteFloatToArray(3, this->_internal_temporal_window(), target);
  }

  // optional string classifiers_property_file_path = 4 [default = "modules/perception/model/sequence_type_fuser/classifiers.property"];
  if (cached_has_bits & 0x00000004u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_classifiers_property_file_path().data(), static_cast<int>(this->_internal_classifiers_property_file_path().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "apollo.perception.sequence_type_fuser_config.ModelConfigs.classifiers_property_file_path");
    target = stream->WriteStringMaybeAliased(
        4, this->_internal_classifiers_property_file_path(), target);
  }

  // optional string transition_property_file_path = 5 [default = "modules/perception/model/sequence_type_fuser/transition.property"];
  if (cached_has_bits & 0x00000008u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_transition_property_file_path().data(), static_cast<int>(this->_internal_transition_property_file_path().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "apollo.perception.sequence_type_fuser_config.ModelConfigs.transition_property_file_path");
    target = stream->WriteStringMaybeAliased(
        5, this->_internal_transition_property_file_path(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:apollo.perception.sequence_type_fuser_config.ModelConfigs)
  return target;
}

size_t ModelConfigs::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:apollo.perception.sequence_type_fuser_config.ModelConfigs)
  size_t total_size = 0;

  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _impl_._has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    // optional string name = 1 [default = "SequenceTypeFuser"];
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
          this->_internal_name());
    }

    // optional string version = 2 [default = "1.1.0"];
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
          this->_internal_version());
    }

    // optional string classifiers_property_file_path = 4 [default = "modules/perception/model/sequence_type_fuser/classifiers.property"];
    if (cached_has_bits & 0x00000004u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
          this->_internal_classifiers_property_file_path());
    }

    // optional string transition_property_file_path = 5 [default = "modules/perception/model/sequence_type_fuser/transition.property"];
    if (cached_has_bits & 0x00000008u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
          this->_internal_transition_property_file_path());
    }

    // optional float temporal_window = 3 [default = 20];
    if (cached_has_bits & 0x00000010u) {
      total_size += 1 + 4;
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
// @@protoc_insertion_point(class_specific_merge_from_start:apollo.perception.sequence_type_fuser_config.ModelConfigs)
  GOOGLE_DCHECK_NE(&from, this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._impl_._has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    if (cached_has_bits & 0x00000001u) {
      _internal_set_name(from._internal_name());
    }
    if (cached_has_bits & 0x00000002u) {
      _internal_set_version(from._internal_version());
    }
    if (cached_has_bits & 0x00000004u) {
      _internal_set_classifiers_property_file_path(from._internal_classifiers_property_file_path());
    }
    if (cached_has_bits & 0x00000008u) {
      _internal_set_transition_property_file_path(from._internal_transition_property_file_path());
    }
    if (cached_has_bits & 0x00000010u) {
      _impl_.temporal_window_ = from._impl_.temporal_window_;
    }
    _impl_._has_bits_[0] |= cached_has_bits;
  }
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void ModelConfigs::CopyFrom(const ModelConfigs& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:apollo.perception.sequence_type_fuser_config.ModelConfigs)
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
      &_impl_.classifiers_property_file_path_, lhs_arena,
      &other->_impl_.classifiers_property_file_path_, rhs_arena
  );
  ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr::InternalSwap(
      &_impl_.transition_property_file_path_, lhs_arena,
      &other->_impl_.transition_property_file_path_, rhs_arena
  );
  swap(_impl_.temporal_window_, other->_impl_.temporal_window_);
}

::PROTOBUF_NAMESPACE_ID::Metadata ModelConfigs::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto_getter, &descriptor_table_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto_once,
      file_level_metadata_modules_2fperception_2fproto_2fsequence_5ftype_5ffuser_5fconfig_2eproto[0]);
}

// @@protoc_insertion_point(namespace_scope)
}  // namespace sequence_type_fuser_config
}  // namespace perception
}  // namespace apollo
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::apollo::perception::sequence_type_fuser_config::ModelConfigs*
Arena::CreateMaybeMessage< ::apollo::perception::sequence_type_fuser_config::ModelConfigs >(Arena* arena) {
  return Arena::CreateMessageInternal< ::apollo::perception::sequence_type_fuser_config::ModelConfigs >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
