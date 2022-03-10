/***************************************************************************
Copyright 2018 The Apollo Authors. All Rights Reserved                     /
                                                                            /
Licensed under the Apache License, Version 2.0 (the "License");             /
you may not use this file except in compliance with the License.            /
You may obtain a copy of the License at                                     /
                                                                            /
    http://www.apache.org/licenses/LICENSE-2.0                              /
                                                                            /
Unless required by applicable law or agreed to in writing, software         /
distributed under the License is distributed on an "AS IS" BASIS,           /
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    /
See the License for the specific language governing permissions and         /
limitations under the License.                                              /
****************************************************************************/

#ifndef LSLIDAR_DECODER_NODELET_H
#define LSLIDAR_DECODER_NODELET_H

#include <nodelet/nodelet.h>
#include <pluginlib/class_list_macros.h>
#include <ros/ros.h>

#include <lslidar_decoder/lslidar_decoder.h>

namespace apollo {
namespace drivers {
namespace lslidar_decoder {
class LslidarDecoderNodelet : public nodelet::Nodelet {
 public:
  LslidarDecoderNodelet() {}
  ~LslidarDecoderNodelet() {}

 private:
  virtual void      onInit();
  LslidarDecoderPtr decoder;
};

}  // namespace lslidar_decoder
}  // namespace drivers
}  // namespace apollo

#endif
