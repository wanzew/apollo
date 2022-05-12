/******************************************************************************
 * Copyright 2018 The Apollo Authors. All Rights Reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *****************************************************************************/

#include "modules/canbus/vehicle/lexus/protocol/wiper_rpt_234.h"

#include "glog/logging.h"

#include "modules/drivers/canbus/common/byte.h"
#include "modules/drivers/canbus/common/canbus_consts.h"

namespace apollo {
namespace canbus {
namespace lexus {

using ::apollo::drivers::canbus::Byte;

Wiperrpt234::Wiperrpt234() {}
const int32_t Wiperrpt234::ID = 0x234;

void Wiperrpt234::Parse(const std::uint8_t* bytes, int32_t length, ChassisDetail* chassis) const {
  chassis->mutable_lexus()->mutable_wiper_rpt_234()->set_vehicle_fault(
      vehicle_fault(bytes, length));
  chassis->mutable_lexus()->mutable_wiper_rpt_234()->set_pacmod_fault(pacmod_fault(bytes, length));
  chassis->mutable_lexus()->mutable_wiper_rpt_234()->set_override_active(
      override_active(bytes, length));
  chassis->mutable_lexus()->mutable_wiper_rpt_234()->set_output_reported_fault(
      output_reported_fault(bytes, length));
  chassis->mutable_lexus()->mutable_wiper_rpt_234()->set_input_output_fault(
      input_output_fault(bytes, length));
  chassis->mutable_lexus()->mutable_wiper_rpt_234()->set_enabled(enabled(bytes, length));
  chassis->mutable_lexus()->mutable_wiper_rpt_234()->set_command_output_fault(
      command_output_fault(bytes, length));
  chassis->mutable_lexus()->mutable_wiper_rpt_234()->set_output_value(output_value(bytes, length));
  chassis->mutable_lexus()->mutable_wiper_rpt_234()->set_commanded_value(
      commanded_value(bytes, length));
  chassis->mutable_lexus()->mutable_wiper_rpt_234()->set_manual_input(manual_input(bytes, length));
}

// config detail: {'name': 'vehicle_fault', 'offset': 0.0, 'precision': 1.0,
// 'len': 1, 'is_signed_var': False, 'physical_range': '[0|1]', 'bit': 6,
// 'type': 'bool', 'order': 'motorola', 'physical_unit': ''}
bool Wiperrpt234::vehicle_fault(const std::uint8_t* bytes, int32_t length) const {
  Byte    t0(bytes + 0);
  int32_t x = t0.get_byte(6, 1);

  bool ret = x;
  return ret;
}

// config detail: {'name': 'pacmod_fault', 'offset': 0.0, 'precision': 1.0,
// 'len': 1, 'is_signed_var': False, 'physical_range': '[0|1]', 'bit': 5,
// 'type': 'bool', 'order': 'motorola', 'physical_unit': ''}
bool Wiperrpt234::pacmod_fault(const std::uint8_t* bytes, int32_t length) const {
  Byte    t0(bytes + 0);
  int32_t x = t0.get_byte(5, 1);

  bool ret = x;
  return ret;
}

// config detail: {'name': 'override_active', 'offset': 0.0, 'precision': 1.0,
// 'len': 1, 'is_signed_var': False, 'physical_range': '[0|1]', 'bit': 1,
// 'type': 'bool', 'order': 'motorola', 'physical_unit': ''}
bool Wiperrpt234::override_active(const std::uint8_t* bytes, int32_t length) const {
  Byte    t0(bytes + 0);
  int32_t x = t0.get_byte(1, 1);

  bool ret = x;
  return ret;
}

// config detail: {'name': 'output_reported_fault', 'offset': 0.0,
// 'precision': 1.0, 'len': 1, 'is_signed_var': False, 'physical_range':
// '[0|1]', 'bit': 4, 'type': 'bool', 'order': 'motorola', 'physical_unit': ''}
bool Wiperrpt234::output_reported_fault(const std::uint8_t* bytes, int32_t length) const {
  Byte    t0(bytes + 0);
  int32_t x = t0.get_byte(4, 1);

  bool ret = x;
  return ret;
}

// config detail: {'name': 'input_output_fault', 'offset': 0.0,
// 'precision': 1.0, 'len': 1, 'is_signed_var': False, 'physical_range':
// '[0|1]', 'bit': 3, 'type': 'bool', 'order': 'motorola', 'physical_unit': ''}
bool Wiperrpt234::input_output_fault(const std::uint8_t* bytes, int32_t length) const {
  Byte    t0(bytes + 0);
  int32_t x = t0.get_byte(3, 1);

  bool ret = x;
  return ret;
}

// config detail: {'name': 'enabled', 'offset': 0.0, 'precision': 1.0, 'len': 1,
// 'is_signed_var': False, 'physical_range': '[0|1]', 'bit': 0, 'type': 'bool',
// 'order': 'motorola', 'physical_unit': ''}
bool Wiperrpt234::enabled(const std::uint8_t* bytes, int32_t length) const {
  Byte    t0(bytes + 0);
  int32_t x = t0.get_byte(0, 1);

  bool ret = x;
  return ret;
}

// config detail: {'name': 'command_output_fault', 'offset': 0.0,
// 'precision': 1.0, 'len': 1, 'is_signed_var': False, 'physical_range':
// '[0|1]', 'bit': 2, 'type': 'bool', 'order': 'motorola', 'physical_unit': ''}
bool Wiperrpt234::command_output_fault(const std::uint8_t* bytes, int32_t length) const {
  Byte    t0(bytes + 0);
  int32_t x = t0.get_byte(2, 1);

  bool ret = x;
  return ret;
}

// config detail: {'name': 'output_value', 'enum': {0:
// 'OUTPUT_VALUE_WIPERS_OFF', 1: 'OUTPUT_VALUE_INTERMITTENT_1', 2:
// 'OUTPUT_VALUE_INTERMITTENT_2', 3: 'OUTPUT_VALUE_INTERMITTENT_3', 4:
// 'OUTPUT_VALUE_INTERMITTENT_4', 5: 'OUTPUT_VALUE_INTERMITTENT_5', 6:
// 'OUTPUT_VALUE_LOW', 7: 'OUTPUT_VALUE_HIGH'}, 'precision': 1.0, 'len': 8,
// 'is_signed_var': False, 'offset': 0.0, 'physical_range': '[0|7]', 'bit': 31,
// 'type': 'enum', 'order': 'motorola', 'physical_unit': ''}
Wiper_rpt_234::Output_valueType Wiperrpt234::output_value(const std::uint8_t* bytes,
                                                          int32_t             length) const {
  Byte    t0(bytes + 3);
  int32_t x = t0.get_byte(0, 8);

  Wiper_rpt_234::Output_valueType ret = static_cast<Wiper_rpt_234::Output_valueType>(x);
  return ret;
}

// config detail: {'name': 'commanded_value', 'enum': {0:
// 'COMMANDED_VALUE_WIPERS_OFF', 1: 'COMMANDED_VALUE_INTERMITTENT_1', 2:
// 'COMMANDED_VALUE_INTERMITTENT_2', 3: 'COMMANDED_VALUE_INTERMITTENT_3', 4:
// 'COMMANDED_VALUE_INTERMITTENT_4', 5: 'COMMANDED_VALUE_INTERMITTENT_5', 6:
// 'COMMANDED_VALUE_LOW', 7: 'COMMANDED_VALUE_HIGH'}, 'precision': 1.0, 'len':
// 8, 'is_signed_var': False, 'offset': 0.0, 'physical_range': '[0|7]', 'bit':
// 23, 'type': 'enum', 'order': 'motorola', 'physical_unit': ''}
Wiper_rpt_234::Commanded_valueType Wiperrpt234::commanded_value(const std::uint8_t* bytes,
                                                                int32_t             length) const {
  Byte    t0(bytes + 2);
  int32_t x = t0.get_byte(0, 8);

  Wiper_rpt_234::Commanded_valueType ret = static_cast<Wiper_rpt_234::Commanded_valueType>(x);
  return ret;
}

// config detail: {'name': 'manual_input', 'enum': {0:
// 'MANUAL_INPUT_WIPERS_OFF', 1: 'MANUAL_INPUT_INTERMITTENT_1', 2:
// 'MANUAL_INPUT_INTERMITTENT_2', 3: 'MANUAL_INPUT_INTERMITTENT_3', 4:
// 'MANUAL_INPUT_INTERMITTENT_4', 5: 'MANUAL_INPUT_INTERMITTENT_5', 6:
// 'MANUAL_INPUT_LOW', 7: 'MANUAL_INPUT_HIGH'}, 'precision': 1.0, 'len': 8,
// 'is_signed_var': False, 'offset': 0.0, 'physical_range': '[0|7]', 'bit': 15,
// 'type': 'enum', 'order': 'motorola', 'physical_unit': ''}
Wiper_rpt_234::Manual_inputType Wiperrpt234::manual_input(const std::uint8_t* bytes,
                                                          int32_t             length) const {
  Byte    t0(bytes + 1);
  int32_t x = t0.get_byte(0, 8);

  Wiper_rpt_234::Manual_inputType ret = static_cast<Wiper_rpt_234::Manual_inputType>(x);
  return ret;
}
}  // namespace lexus
}  // namespace canbus
}  // namespace apollo
