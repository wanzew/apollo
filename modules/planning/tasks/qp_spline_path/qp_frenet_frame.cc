/******************************************************************************
 * Copyright 2017 The Apollo Authors. All Rights Reserved.
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

/**
 * @file qp_frenet_frame.cc
 **/
#include "modules/planning/tasks/qp_spline_path/qp_frenet_frame.h"

#include <algorithm>
#include <iterator>
#include <limits>

#include "modules/common/proto/pnc_point.pb.h"
#include "modules/planning/proto/planning.pb.h"

#include "modules/common/configs/vehicle_config_helper.h"
#include "modules/common/macro.h"
#include "modules/common/math/linear_interpolation.h"
#include "modules/common/util/util.h"
#include "modules/planning/common/planning_gflags.h"
#include "modules/planning/common/planning_util.h"

namespace apollo {
namespace planning {

using common::SpeedPoint;

namespace {

constexpr double kEpsilonTol = 1e-6;
const auto       inf         = std::numeric_limits<double>::infinity();
}  // namespace

QpFrenetFrame::QpFrenetFrame(const ReferenceLine&            reference_line,
                             const SpeedData&                speed_data,
                             const common::FrenetFramePoint& init_frenet_point,
                             const double                    time_resolution,
                             const std::vector<double>&      evaluated_s)
    : reference_line_(reference_line)
    , speed_data_(speed_data)
    , vehicle_param_(common::VehicleConfigHelper::instance()->GetConfig().vehicle_param())
    , init_frenet_point_(init_frenet_point)
    , feasible_longitudinal_upper_bound_(reference_line_.map_path().length())
    , start_s_(init_frenet_point.s())
    , end_s_(reference_line.Length())
    , time_resolution_(time_resolution)
    , evaluated_s_(evaluated_s)
    , hdmap_bound_{evaluated_s_.size(), std::make_pair(-inf, inf)}
    , static_obstacle_bound_{evaluated_s_.size(), std::make_pair(-inf, inf)}
    , dynamic_obstacle_bound_{evaluated_s_.size(), std::make_pair(-inf, inf)} {
  DCHECK_LT(start_s_ + kEpsilonTol, end_s_);
  DCHECK_GE(evaluated_s_.size(), 2);
}

bool QpFrenetFrame::Init(const std::vector<const PathObstacle*>& path_obstacles) {
  // 根据SpeedData，计算自车纵向行驶轨迹，存入vector<SpeedPoint>
  // 即discretized_vehicle_location_
  // calculate discretized vehicle location
  CalculateDiscretizedVehicleLocation();

  //根据HDMap和道路参考线，计算自车轨迹的横向约束
  // Calculate HDMap bound
  CalculateHDMapBound();

  // 计算静态、动态障碍物及相应的decision对自车轨迹的横向约束
  // Calculate obstacle bound
  CalculateObstacleBound(path_obstacles);

  return true;
}

void QpFrenetFrame::LogQpBound(apollo::planning_internal::Debug* planning_debug) {
  if (!planning_debug) { return; }
  apollo::planning_internal::SLFrameDebug* sl_frame =
      planning_debug->mutable_planning_data()->mutable_sl_frame()->Add();
  for (size_t i = 0; i < evaluated_s_.size(); ++i) {
    sl_frame->mutable_sampled_s()->Add(evaluated_s_[i]);
    sl_frame->mutable_map_lower_bound()->Add(hdmap_bound_[i].first);
    sl_frame->mutable_map_upper_bound()->Add(hdmap_bound_[i].second);
    sl_frame->mutable_static_obstacle_lower_bound()->Add(static_obstacle_bound_[i].first);
    sl_frame->mutable_static_obstacle_upper_bound()->Add(static_obstacle_bound_[i].second);
    sl_frame->mutable_dynamic_obstacle_lower_bound()->Add(dynamic_obstacle_bound_[i].first);
    sl_frame->mutable_dynamic_obstacle_upper_bound()->Add(dynamic_obstacle_bound_[i].second);
  }
}

const std::vector<std::pair<double, double>>& QpFrenetFrame::GetMapBound() const {
  return hdmap_bound_;
}

const std::vector<std::pair<double, double>>& QpFrenetFrame::GetStaticObstacleBound() const {
  return static_obstacle_bound_;
}

const std::vector<std::pair<double, double>>& QpFrenetFrame::GetDynamicObstacleBound() const {
  return dynamic_obstacle_bound_;
}

// 根据SpeedData，计算自车纵向行经轨迹，存入vector<SpeedPoint>
// 即discretized_vehicle_location_
bool QpFrenetFrame::CalculateDiscretizedVehicleLocation() {
  // 根据SpeedData按时间遍历计算了自车的纵向轨迹相关信息，
  // 存入vector<SpeedPoint>，即discretized_vehicle_location_。
  for (double relative_time = 0.0; relative_time < speed_data_.TotalTime();
       relative_time += time_resolution_) {
    SpeedPoint veh_point;
    // 此处的speed_data_是从 DP_SPEED 计算得到的结果，并不一定是匀速的了
    speed_data_.EvaluateByTime(relative_time, &veh_point);
    veh_point.set_t(relative_time);
    discretized_vehicle_location_.push_back(std::move(veh_point));
  }
  return true;
}

// MapDynamicObstacleWithDecision()也是和静态障碍物大体相似的思路。
// 根据由SpeedData求得的自车位置信息discretized_vehicle_location_，
// 按照时间戳一一匹配动态障碍物的轨迹点（应该是预测轨迹），计算该点处的障碍物bounding
// box对自车横向轨迹的影响。
bool QpFrenetFrame::MapDynamicObstacleWithDecision(const PathObstacle& path_obstacle) {
  const Obstacle* ptr_obstacle = path_obstacle.obstacle();
  if (!path_obstacle.HasLateralDecision()) {
    ADEBUG << "object has no lateral decision";
    return false;
  }
  const auto& decision = path_obstacle.LateralDecision();
  if (!decision.has_nudge()) {
    ADEBUG << "only support nudge now";
    return true;
  }
  const auto& nudge = decision.nudge();
  for (const SpeedPoint& veh_point : discretized_vehicle_location_) {
    double                  time             = veh_point.t();
    common::TrajectoryPoint trajectory_point = ptr_obstacle->GetPointAtTime(time);
    common::math::Box2d     obs_box          = ptr_obstacle->GetBoundingBox(trajectory_point);
    // project obs_box on reference line
    std::vector<common::math::Vec2d> corners;
    obs_box.GetAllCorners(&corners);
    std::vector<common::SLPoint> sl_corners;

    for (const auto& corner_xy : corners) {
      common::SLPoint cur_point;
      reference_line_.XYToSL(corner_xy, &cur_point);
      // shift box base on buffer
      cur_point.set_l(cur_point.l() + nudge.distance_l());
      sl_corners.push_back(std::move(cur_point));
    }

    for (uint32_t i = 0; i < sl_corners.size(); ++i) {
      common::SLPoint sl_first  = sl_corners[i % sl_corners.size()];
      common::SLPoint sl_second = sl_corners[(i + 1) % sl_corners.size()];
      if (sl_first.s() < sl_second.s()) { std::swap(sl_first, sl_second); }

      // MapLateralConstraint()被调用4次，因为sl_corners由bounding box而来，有4个端点
      //因为是考虑横向约束，bounding box的上下2条横向线段端点对于求bound没有影响
      std::pair<double, double> bound = MapLateralConstraint(      //
          sl_first,                                                //
          sl_second,                                               //
          nudge.type(),                                            //
          veh_point.s() - vehicle_param_.back_edge_to_center(),    //
          veh_point.s() + vehicle_param_.front_edge_to_center());  //

      // update bound map
      double s_resolution    = std::fabs(veh_point.v() * time_resolution_);
      double updated_start_s = init_frenet_point_.s() + veh_point.s() - s_resolution;
      double updated_end_s   = init_frenet_point_.s() + veh_point.s() + s_resolution;
      // 纵向s超出考察范围
      if (updated_end_s > evaluated_s_.back() || updated_start_s < evaluated_s_.front()) {
        continue;
      }
      // clang-format off
      std::pair<uint32_t, uint32_t> update_index_range = FindInterval(updated_start_s, updated_end_s);
      for (uint32_t j = update_index_range.first; j <= update_index_range.second; ++j) {
        //可行驶区域的右边界bound，取max
        dynamic_obstacle_bound_[j].first = std::max(bound.first, dynamic_obstacle_bound_[j].first);
        //可行驶区域的左边界bound，取min
        dynamic_obstacle_bound_[j].second = std::min(bound.second, dynamic_obstacle_bound_[j].second);
      }
      // clang-format on
    }
  }
  return true;
}

// 考虑了静态障碍物以及相应的横向避让措施对可行驶区域宽度的影响，
// 处理结果static_obstacle_bound_保存了evaluated_s_中各采样点处的横向可行驶范围。
bool QpFrenetFrame::MapStaticObstacleWithDecision(const PathObstacle& path_obstacle) {
  const auto ptr_obstacle = path_obstacle.obstacle();
  if (!path_obstacle.HasLateralDecision()) {
    ADEBUG << "obstacle has no lateral decision";
    return false;
  }
  const auto& decision = path_obstacle.LateralDecision();
  if (!decision.has_nudge()) {
    ADEBUG << "only support nudge decision now";
    return true;
  }
  // 处理静态障碍物的主要思路是将其轮廓端点映射到Frenet坐标系，
  // 结合nudge的方向计算障碍物占据的横向范围，由 MapNudgePolygon()和 MapNudgeLine()实现。
  if (!MapNudgePolygon(common::math::Polygon2d(ptr_obstacle->PerceptionBoundingBox()),
                       decision.nudge(), &static_obstacle_bound_)) {
    AERROR << "fail to map polygon with id " << path_obstacle.Id() << " in qp frenet frame";
    return false;
  }
  return true;
}

//确定障碍物轮廓对纵向s范围内轨迹采样点的横向约束
bool QpFrenetFrame::MapNudgePolygon(const common::math::Polygon2d&                polygon,
                                    const ObjectNudge&                            nudge,
                                    std::vector<std::pair<double, double>>* const bound_map) {
  std::vector<common::SLPoint> sl_corners;
  for (const auto& corner_xy : polygon.points()) {
    common::SLPoint corner_sl;
    if (!reference_line_.XYToSL(corner_xy, &corner_sl)) {
      AERROR << "Fail to map xy point " << corner_xy.DebugString() << " to "
             << corner_sl.DebugString();
      return false;
    }
    // shift box based on buffer
    // nudge decision buffer:
    // --- position for left nudge
    // --- negative for right nudge
    corner_sl.set_l(corner_sl.l() + nudge.distance_l());
    sl_corners.push_back(std::move(corner_sl));
  }

  const auto corner_size = sl_corners.size();
  // MapNudgePolygon()中循环调用MapNudgeLine()，
  // 将障碍物的bounding box对自车轨迹的横向约束计算，转换为计算bounding
  // box的4条边对自车轨迹的横向约束。
  for (uint32_t i = 0; i < corner_size; ++i) {
    if (!MapNudgeLine(sl_corners[i], sl_corners[(i + 1) % corner_size], nudge.type(), bound_map)) {
      AERROR << "Map box line (sl) " << sl_corners[i].DebugString() << "->"
             << sl_corners[(i + 1) % corner_size].DebugString();
      return false;
    }
  }
  return true;
}

//确定障碍物轮廓的一条边所对应纵向s范围内轨迹采样点的横向约束
bool QpFrenetFrame::MapNudgeLine(const common::SLPoint&                        start,
                                 const common::SLPoint&                        end,
                                 const ObjectNudge::Type                       nudge_type,
                                 std::vector<std::pair<double, double>>* const constraint) {
  DCHECK_NOTNULL(constraint);

  const common::SLPoint& near_point    = (start.s() < end.s() ? start : end);
  const common::SLPoint& further_point = (start.s() < end.s() ? end : start);
  // impact_index表示evaluated_s_中受影响的区间范围，2个index可能是相等的
  std::pair<uint32_t, uint32_t> impact_index = FindInterval(near_point.s(), further_point.s());
  // s轴超出轨迹纵向范围
  if (further_point.s() < start_s_ - vehicle_param_.back_edge_to_center() ||
      near_point.s() > end_s_ + vehicle_param_.front_edge_to_center()) {
    return true;
  }

  const double distance = std::max(further_point.s() - near_point.s(), common::math::kMathEpsilon);
  const double adc_half_width = vehicle_param_.width() / 2;

  DCHECK_GT(constraint->size(), impact_index.second);
  for (uint32_t i = impact_index.first; i <= impact_index.second; ++i) {
    double weight   = std::fabs((evaluated_s_[i] - near_point.s())) / distance;
    weight          = std::min(1.0, std::max(weight, 0.0));
    double boundary = near_point.l() * (1 - weight) + further_point.l() * weight;

    if (nudge_type == ObjectNudge::LEFT_NUDGE) {
      boundary += adc_half_width;
      // first是lower bound，second是upper bound
      // lower bound增大，即将自车可行驶区域向左移动，即left nudge
      //而upper bound不变，其初始化为INF
      (*constraint)[i].first = std::max(boundary, (*constraint)[i].first);
    } else {
      boundary -= adc_half_width;
      (*constraint)[i].second = std::min(boundary, (*constraint)[i].second);
    }
    //若可行驶区域宽度constraint太窄(<0.3m)，则收缩feasible_longitudinal_upper_bound_
    if ((*constraint)[i].second < (*constraint)[i].first + 0.3) {
      if (i > 0) {
        feasible_longitudinal_upper_bound_ =
            std::min(evaluated_s_[i - 1] - kEpsilonTol, feasible_longitudinal_upper_bound_);
      } else {
        feasible_longitudinal_upper_bound_ = start_s_;
        return true;
      }

      ADEBUG << "current mapping constraint, sl point impact index "
             << "near_point: " << near_point.DebugString()
             << "further_point: " << further_point.DebugString() << "impact_index: " << impact_index
             << "(*constraint)[" << i << "]" << (*constraint)[i];
      break;
    }
  }

  return true;
}

// 用来计算障碍物bounding box的一条边对自车横向轨迹的约束。
//返回值pair.first表示可行驶区域的右边界限定，second表示左边界限定，second > first
std::pair<double, double> QpFrenetFrame::MapLateralConstraint(  //
    const common::SLPoint&  start,
    const common::SLPoint&  end,
    const ObjectNudge::Type nudge_type,
    const double            s_start,
    const double            s_end) {
  constexpr double          inf    = std::numeric_limits<double>::infinity();
  std::pair<double, double> result = std::make_pair(-inf, inf);

  //障碍物车在自车前方或后方，忽略，本函数只考虑横向影响
  if (start.s() > s_end || end.s() < s_start) { return result; }
  double s_front = std::max(start.s(), s_start);
  double s_back  = std::min(end.s(), s_end);

  double weight_back  = 0.0;
  double weight_front = 0.0;

  if (end.s() - start.s() > std::numeric_limits<double>::epsilon()) {
    weight_back  = (s_back - end.s()) / (end.s() - start.s());
    weight_front = (s_front - start.s()) / (end.s() - start.s());
  }

  //将自车首尾点向障碍物车做映射，找对应点，以确定横向偏移
  //我觉得只简单的利用障碍物车的bounding box就够了
  // nudge的横向距离已经在调用该函数之前纳入障碍物车的start和end两点的l坐标了
  // clang-format off
  common::SLPoint front =common::math::InterpolateUsingLinearApproximation(start, end, weight_front);
  common::SLPoint back = common::math::InterpolateUsingLinearApproximation(start, end, weight_back);
  // clang-format on

  if (nudge_type == ObjectNudge::RIGHT_NUDGE) {
    //确认自车横向左方限定，取min使可行驶范围偏右，故nudge right
    result.second = std::min(front.l(), back.l());
  } else {
    //确认自车横向右方限定，取max使可行驶范围偏左，故nudge left
    result.first = std::max(front.l(), back.l());
  }
  return result;
}

std::pair<uint32_t, uint32_t> QpFrenetFrame::FindInterval(const double start,
                                                          const double end) const {
  double new_start = std::max(start - vehicle_param_.front_edge_to_center(), evaluated_s_.front());
  double new_end   = std::min(end + vehicle_param_.back_edge_to_center(), evaluated_s_.back());

  uint32_t start_index = FindIndex(new_start);
  uint32_t end_index   = FindIndex(new_end);
  if (end_index == evaluated_s_.size()) { --end_index; }

  return std::make_pair(start_index, end_index);
}

// 计算了沿纵轴s均匀采样所得点集中 各个点的横向约束，即左右边界。
bool QpFrenetFrame::CalculateHDMapBound() {
  // 根据HDMap和道路参考线，计算自车轨迹的横向约束
  // hdmap_bound_初始化为vector {evaluated_s_.size(), std::make_pair(-inf, inf)}
  const double adc_half_width = vehicle_param_.width() / 2.0;
  for (uint32_t i = 0; i < hdmap_bound_.size(); ++i) {
    double left_bound  = 0.0;
    double right_bound = 0.0;
    bool   suc         = reference_line_.GetLaneWidth(evaluated_s_[i], &left_bound, &right_bound);
    if (!suc) {
      AWARN << "Extracting lane width failed at s = " << evaluated_s_[i];
      right_bound = FLAGS_default_reference_line_width / 2;
      left_bound  = FLAGS_default_reference_line_width / 2;
    }

    // 按照右手坐标系，自车前方为s轴正方向，
    // hdmap_bound_[i].first是右侧边界，second是左侧边界， 左右各往里收半个车宽
    hdmap_bound_[i].first  = -right_bound + adc_half_width;
    hdmap_bound_[i].second = left_bound - adc_half_width;

    // 如果右边界>左边界，不合理，则缩短纵向s上限
    // feasible_longitudinal_upper_bound_到当前考察点
    if (hdmap_bound_[i].first >= hdmap_bound_[i].second) {
      feasible_longitudinal_upper_bound_ =
          std::min(evaluated_s_[i], feasible_longitudinal_upper_bound_);
      common::SLPoint sl;
      sl.set_s(evaluated_s_[i]);

      common::math::Vec2d xy;
      reference_line_.SLToXY(sl, &xy);
      ADEBUG << "evaluated point x: " << std::fixed << xy.x() << " y: " << xy.y();
      break;
    }
  }
  return true;
}

// CalculateObstacleBound()的处理分为2部分：静态障碍物和动态障碍物。
// 首先判断障碍物是否有LateralDecision，若没有，则忽略，因为对path没有影响（这个阶段的path重点考虑横向运动）。
// LateralDecision主要是指横向的nudge，
// 具体decision相关信息可查阅apollo3_5/modules/planning/proto/decision.proto。
bool QpFrenetFrame::CalculateObstacleBound(const std::vector<const PathObstacle*>& path_obstacles) {
  for (const auto ptr_path_obstacle : path_obstacles) {
    // path是横向轨迹，若没有LateralDecision，就对path优化没有影响，主要指nudge
    if (!ptr_path_obstacle->HasLateralDecision()) { continue; }
    if (ptr_path_obstacle->obstacle()->IsStatic()) {
      //计算静态、动态障碍物及相应的decision对自车轨迹的横向约束
      if (!MapStaticObstacleWithDecision(*ptr_path_obstacle)) {
        AERROR << "mapping obstacle with id [" << ptr_path_obstacle->Id()
               << "] failed in qp frenet frame.";
        return false;
      }
    } else {
      if (!MapDynamicObstacleWithDecision(*ptr_path_obstacle)) {
        AERROR << "mapping obstacle with id [" << ptr_path_obstacle->Id()
               << "] failed in qp frenet frame.";
        return false;
      }
    }
  }
  return true;
}

uint32_t QpFrenetFrame::FindIndex(const double s) const {
  auto upper_bound = std::upper_bound(evaluated_s_.begin(), evaluated_s_.end(), s);
  return std::distance(evaluated_s_.begin(), upper_bound);
}

}  // namespace planning
}  // namespace apollo
