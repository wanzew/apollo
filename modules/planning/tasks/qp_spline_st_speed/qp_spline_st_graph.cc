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
 * @file
 **/

#include "modules/planning/tasks/qp_spline_st_speed/qp_spline_st_graph.h"

#include <algorithm>
#include <limits>
#include <string>
#include <utility>

#include "modules/common/log.h"
#include "modules/planning/common/frame.h"
#include "modules/planning/common/planning_gflags.h"

namespace apollo {
namespace planning {

using apollo::common::ErrorCode;
using apollo::common::SpeedPoint;
using apollo::common::Status;
using apollo::common::VehicleParam;
using apollo::planning_internal::STGraphDebug;

QpSplineStGraph::QpSplineStGraph(Spline1dGenerator*     spline_generator,
                                 const QpStSpeedConfig& qp_st_speed_config,
                                 const VehicleParam&    veh_param,
                                 const bool             is_change_lane)
    : spline_generator_(spline_generator)
    , qp_st_speed_config_(qp_st_speed_config)
    , is_change_lane_(is_change_lane)
    , t_knots_resolution_(qp_st_speed_config_.total_time() /  // 默认时间长度 8.0 s
                          qp_st_speed_config_.qp_spline_config().number_of_discrete_graph_t())  // 4
{
  Init();
}

void QpSplineStGraph::Init() {
  // init knots
  double   curr_t     = 0.0;
  uint32_t num_spline = qp_st_speed_config_.qp_spline_config().number_of_discrete_graph_t() - 1;
  for (uint32_t i = 0; i <= num_spline; ++i) {
    t_knots_.push_back(curr_t);
    curr_t += t_knots_resolution_;
  }

  uint32_t num_evaluated_t = 10 * num_spline + 1;
  curr_t                   = 0;  // init evaluated t positions
  t_evaluated_resolution_  = qp_st_speed_config_.total_time() / (num_evaluated_t - 1);

  for (uint32_t i = 0; i < num_evaluated_t; ++i) {
    t_evaluated_.push_back(curr_t);
    curr_t += t_evaluated_resolution_;
  }
}

void QpSplineStGraph::SetDebugLogger(planning_internal::STGraphDebug* st_graph_debug) {
  if (st_graph_debug) {
    st_graph_debug->Clear();
    st_graph_debug_ = st_graph_debug;
  }
}

// 规划的主体过程在QpSplineStGraph::Search()中，分为4个步骤：
// AddConstraint()，
// AddKernel()，
// Solve()，
// extract speed data。
Status QpSplineStGraph::Search(const StGraphData&               st_graph_data,
                               const std::pair<double, double>& accel_bound,
                               const SpeedData&                 reference_speed_data,
                               SpeedData* const                 speed_data) {
  constexpr double kBounadryEpsilon = 1e-2;
  for (auto boundary : st_graph_data.st_boundaries()) {
    if (boundary->IsPointInBoundary({0.0, 0.0}) ||
        (std::fabs(boundary->min_t()) < kBounadryEpsilon &&
         std::fabs(boundary->min_s()) < kBounadryEpsilon)) {
      speed_data->Clear();
      const double t_output_resolution = FLAGS_trajectory_time_min_interval;
      double       time                = 0.0;
      while (time < qp_st_speed_config_.total_time() + t_output_resolution) {
        speed_data->AppendSpeedPoint(0.0, time, 0.0, 0.0, 0.0);
        time += t_output_resolution;
      }
      return Status::OK();
    }
  }

  cruise_.clear();
  reference_dp_speed_points_ = speed_data->speed_vector();  // DP 输出的结果
  init_point_                = st_graph_data.init_point() spline_generator_->Reset(
      t_knots_, qp_st_speed_config_.qp_spline_config().spline_order());

  AddConstraint(st_graph_data.init_point(), st_graph_data.speed_limit(),
                st_graph_data.st_boundaries(), accel_bound);
  AddKernel(st_graph_data.st_boundaries(), st_graph_data.speed_limit());
  Solve();
  // extract output
  speed_data->Clear();
  const Spline1d& spline = spline_generator_->spline();

  const double t_output_resolution = FLAGS_trajectory_time_min_interval;
  double       time                = 0.0;
  while (time < qp_st_speed_config_.total_time() + t_output_resolution) {
    double s  = spline(time);
    double v  = std::max(0.0, spline.Derivative(time));
    double a  = spline.SecondOrderDerivative(time);
    double da = spline.ThirdOrderDerivative(time);
    speed_data->AppendSpeedPoint(s, time, v, a, da);
    time += t_output_resolution;
  }

  return Status::OK();
}

Status QpSplineStGraph::AddConstraint(const common::TrajectoryPoint&        init_point,
                                      const SpeedLimit&                     speed_limit,
                                      const std::vector<const StBoundary*>& boundaries,
                                      const std::pair<double, double>&      accel_bound) {
  // clang-format off
  // 主要添加了如下约束
  // 1. 起始点函数值约束 f(0) = 0。等式约束。
  // 2. 起始点速度约束 f'(0) = v0 = 初始速度。等式约束。
  // 3. 单调性约束，即随时间t 增大，s = f(t) 一定是不变或增大（单调递增），因为车会停下或前进，不后退。不等式约束。
  // 4. 分界点处各阶导数连续平滑性约束，这里用了0~3阶导数。等式约束。
  // 5. 采样时刻t 对应的自车s 坐标范围约束。不等式约束。
  // 6. 采样时刻t 对应的自车速度范围约束。不等式约束。
  // 7. 采样时刻t 对应的自车加速度范围约束。不等式约束。
  // clang-format on
  Spline1dConstraint* constraint = spline_generator_->mutable_spline_constraint();

  // add st start point constraint
  constraint->AddPointConstraint(0.0, 0.0);
  // add st start point velocity constrain
  constraint->AddPointDerivativeConstraint(0.0, init_point_.v());
  // add monotone inequality constraint
  constraint->AddMonotoneInequalityConstraint(t_evaluated_);
  // add smoothness joint constraint
  constraint->AddThirdDerivativeSmoothConstraint();

  // boundary constraint
  std::vector<double> s_upper_bound, s_lower_bound;

  // clang-format off
  //求不同时刻t，自车s应该处的范围，每一个t对应一个[lower_s, upper_s]
  for (const double curr_t : t_evaluated_) {
    double lower_s = 0.0;
    double upper_s = 0.0;
    // 根据在time时刻各个boundary的类别，确定s的上界或下界
    // 对应时刻 curr_t 所对应的自车可行驶区域
    // s_lower_bound <= boundaries <= s_upper_bound
    GetSConstraintByTime(boundaries, 
                         curr_t, 
                         qp_st_speed_config_.total_path_length(), 
                         &upper_s,
                         &lower_s);
    s_upper_bound.push_back(upper_s);
    s_lower_bound.push_back(lower_s);
  }
  // clang-format on
  DCHECK_EQ(t_evaluated_.size(), s_lower_bound.size());
  DCHECK_EQ(t_evaluated_.size(), s_upper_bound.size());

  // apply distance constraints
  //添加lower_s <= f(t) <= upper_s取值范围约束
  constraint->AddBoundary(t_evaluated_, s_lower_bound, s_upper_bound);

  /*
   *     |s/m
   *     |
   * 149 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   * 148 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   * 147 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   * 146 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   * 145 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   * 144 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   . |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   . |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   . |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   . |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   . |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   . |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   . |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   . |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *  20 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *  19 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *  18 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *  17 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *  16 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *  15 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *  14 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *  13 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *  12 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *  11 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *  10 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   9 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   8 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   6 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   5 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   4 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   3 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   2 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   1 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *   0 |  (0，s) (1，s) (2，s) (3，s) (4，s) (5，s) (6，s) (7，s)
   *     |________________________________________________________
   *          0      1      2      3     4     5      6      7  t/s
   */
  // speed constraint
  std::vector<double> speed_upper_bound;
  // estimate speed upper constraints
  EstimateSpeedUpperBound(init_point, speed_limit, &speed_upper_bound);

  std::vector<double> speed_lower_bound(t_evaluated_.size(), 0.0);
  DCHECK_EQ(t_evaluated_.size(), speed_upper_bound.size());
  DCHECK_EQ(t_evaluated_.size(), speed_lower_bound.size());

  if (st_graph_debug_) {
    auto speed_constraint = st_graph_debug_->mutable_speed_constraint();
    for (size_t i = 0; i < t_evaluated_.size(); ++i) {
      speed_constraint->add_t(t_evaluated_[i]);
      speed_constraint->add_lower_bound(speed_lower_bound[i]);
      speed_constraint->add_upper_bound(speed_upper_bound[i]);
    }
  }
  // apply speed constraints
  constraint->AddDerivativeBoundary(t_evaluated_, speed_lower_bound, speed_upper_bound);

  // acceleration constraint
  std::vector<double> accel_lower_bound(t_evaluated_.size(), accel_bound.first);
  std::vector<double> accel_upper_bound(t_evaluated_.size(), accel_bound.second);

  bool   has_follow = false;
  double delta_s    = 1.0;
  for (const auto* boundary : boundaries) {
    if (boundary->boundary_type() == StBoundary::BoundaryType::FOLLOW) {
      has_follow = true;
      delta_s    = std::fmin(delta_s, boundary->min_s() - fabs(boundary->characteristic_length()));
    }
  }

  //如果在跟车（即前方有障碍物）且距前车很近，设定第1个时间间隔内加速度上限为0
  if (FLAGS_enable_follow_accel_constraint && has_follow && delta_s < 0.0) {
    accel_upper_bound.front() = 0.0;
  }

  DCHECK_EQ(t_evaluated_.size(), accel_lower_bound.size());
  DCHECK_EQ(t_evaluated_.size(), accel_upper_bound.size());
  // apply acceleration constraints
  //添加accel_lower_bound <= f''(t) <= accel_upper_bound取值范围约束
  constraint->AddSecondDerivativeBoundary(t_evaluated_, accel_lower_bound, accel_upper_bound);

  return Status::OK();
}

Status QpSplineStGraph::AddKernel(const std::vector<const StBoundary*>& boundaries,
                                  const SpeedLimit&                     speed_limit) {
  Spline1dKernel* spline_kernel = spline_generator_->mutable_spline_kernel();

  if (qp_st_speed_config_.qp_spline_config().accel_kernel_weight() > 0) {
    spline_kernel->AddSecondOrderDerivativeMatrix(
        qp_st_speed_config_.qp_spline_config().accel_kernel_weight());
  }

  if (qp_st_speed_config_.qp_spline_config().jerk_kernel_weight() > 0) {
    spline_kernel->AddThirdOrderDerivativeMatrix(
        qp_st_speed_config_.qp_spline_config().jerk_kernel_weight());
  }

  AddCruiseReferenceLineKernel(qp_st_speed_config_.qp_spline_config().cruise_weight();
  AddFollowReferenceLineKernel(boundaries,
                                    qp_st_speed_config_.qp_spline_config().follow_weight());
  AddYieldReferenceLineKernel(boundaries,
                                   qp_st_speed_config_.qp_spline_config().yield_weight());
  AddDpStReferenceKernel(qp_st_speed_config_.qp_spline_config().dp_st_reference_weight());

  // init point jerk continuous kernel
  (*spline_kernel->mutable_kernel_matrix())(2, 2) +=
      2.0 * 4.0 * qp_st_speed_config_.qp_spline_config().init_jerk_kernel_weight();
  (*spline_kernel->mutable_offset())(2, 0) +=
      -4.0 * init_point_.a() * qp_st_speed_config_.qp_spline_config().init_jerk_kernel_weight();

  spline_kernel->AddRegularization(qp_st_speed_config_.qp_spline_config().regularization_weight());

  return Status::OK();
}

Status QpSplineStGraph::Solve() {
  return spline_generator_->Solve() ? Status::OK() :
                                      Status(ErrorCode::PLANNING_ERROR, "QpSplineStGraph::solve");
}

Status QpSplineStGraph::AddCruiseReferenceLineKernel(const double weight) {
  auto*  spline_kernel = spline_generator_->mutable_spline_kernel();
  double dist_ref      = qp_st_speed_config_.total_path_length();
  for (uint32_t i = 0; i < t_evaluated_.size(); ++i) {
    cruise_.push_back(dist_ref);
  }
  if (st_graph_debug_) {
    auto kernel_cruise_ref = st_graph_debug_->mutable_kernel_cruise_ref();
    kernel_cruise_ref->mutable_t()->Add(t_evaluated_[0]);
    kernel_cruise_ref->mutable_cruise_line_s()->Add(dist_ref);
    for (uint32_t i = 1; i < t_evaluated_.size(); ++i) {
      kernel_cruise_ref->mutable_t()->Add(t_evaluated_[i]);
      kernel_cruise_ref->mutable_cruise_line_s()->Add(cruise_[i]);
    }
  }
  DCHECK_EQ(t_evaluated_.size(), cruise_.size());

  for (std::size_t i = 0; i < t_evaluated_.size(); ++i) {
    ADEBUG << "Cruise Ref S: " << cruise_[i] << " Relative time: " << t_evaluated_[i] << std::endl;
  }

  if (t_evaluated_.size() > 0) {
    spline_kernel->AddReferenceLineKernelMatrix(
        t_evaluated_, cruise_, weight * qp_st_speed_config_.total_time() / t_evaluated_.size());
  }

  return Status::OK();
}

Status
QpSplineStGraph::AddFollowReferenceLineKernel(const std::vector<const StBoundary*>& boundaries,
                                              const double                          weight) {
  auto*               spline_kernel = spline_generator_->mutable_spline_kernel();
  std::vector<double> ref_s;
  std::vector<double> filtered_evaluate_t;
  for (size_t i = 0; i < t_evaluated_.size(); ++i) {
    const double curr_t  = t_evaluated_[i];
    double       s_min   = std::numeric_limits<double>::infinity();
    bool         success = false;
    for (const auto* boundary : boundaries) {
      if (boundary->boundary_type() != StBoundary::BoundaryType::FOLLOW) { continue; }
      if (curr_t < boundary->min_t() || curr_t > boundary->max_t()) { continue; }
      double s_upper = 0.0;
      double s_lower = 0.0;
      if (boundary->GetUnblockSRange(curr_t, &s_upper, &s_lower)) {
        success = true;
        s_min   = std::min(s_min, s_upper - boundary->characteristic_length() -
                                    qp_st_speed_config_.qp_spline_config().follow_drag_distance());
      }
    }
    if (success && s_min < cruise_[i]) {
      filtered_evaluate_t.push_back(curr_t);
      ref_s.push_back(s_min);
      if (st_graph_debug_) {
        auto kernel_follow_ref = st_graph_debug_->mutable_kernel_follow_ref();
        kernel_follow_ref->mutable_t()->Add(curr_t);
        kernel_follow_ref->mutable_follow_line_s()->Add(s_min);
      }
    }
  }
  DCHECK_EQ(filtered_evaluate_t.size(), ref_s.size());

  if (!ref_s.empty()) {
    spline_kernel->AddReferenceLineKernelMatrix(filtered_evaluate_t, ref_s,
                                                weight * qp_st_speed_config_.total_time() /
                                                    t_evaluated_.size());
  }

  for (std::size_t i = 0; i < filtered_evaluate_t.size(); ++i) {
    ADEBUG << "Follow Ref S: " << ref_s[i] << " Relative time: " << filtered_evaluate_t[i]
           << std::endl;
  }
  return Status::OK();
}

Status
QpSplineStGraph::AddYieldReferenceLineKernel(const std::vector<const StBoundary*>& boundaries,
                                             const double                          weight) {
  auto*               spline_kernel = spline_generator_->mutable_spline_kernel();
  std::vector<double> ref_s;
  std::vector<double> filtered_evaluate_t;
  for (size_t i = 0; i < t_evaluated_.size(); ++i) {
    const double curr_t  = t_evaluated_[i];
    double       s_min   = std::numeric_limits<double>::infinity();
    bool         success = false;
    for (const auto* boundary : boundaries) {
      if (boundary->boundary_type() != StBoundary::BoundaryType::YIELD) { continue; }
      if (curr_t < boundary->min_t() || curr_t > boundary->max_t()) { continue; }
      double s_upper = 0.0;
      double s_lower = 0.0;
      if (boundary->GetUnblockSRange(curr_t, &s_upper, &s_lower)) {
        success = true;
        s_min   = std::min(s_min, s_upper - boundary->characteristic_length() -
                                    qp_st_speed_config_.qp_spline_config().yield_drag_distance());
      }
    }
    if (success && s_min < cruise_[i]) {
      filtered_evaluate_t.push_back(curr_t);
      ref_s.push_back(s_min);
    }
  }
  DCHECK_EQ(filtered_evaluate_t.size(), ref_s.size());

  if (!ref_s.empty()) {
    spline_kernel->AddReferenceLineKernelMatrix(filtered_evaluate_t, ref_s,
                                                weight * qp_st_speed_config_.total_time() /
                                                    t_evaluated_.size());
  }

  for (std::size_t i = 0; i < filtered_evaluate_t.size(); ++i) {
    ADEBUG << "Yield Ref S: " << ref_s[i] << " Relative time: " << filtered_evaluate_t[i]
           << std::endl;
  }
  return Status::OK();
}

bool QpSplineStGraph::AddDpStReferenceKernel(const double weight) const {
  std::vector<double> t_pos;
  std::vector<double> s_pos;
  for (auto point : reference_dp_speed_points_) {
    t_pos.push_back(point.t());
    s_pos.push_back(point.s());
  }
  auto* spline_kernel = spline_generator_->mutable_spline_kernel();
  if (!t_pos.empty()) {
    spline_kernel->AddReferenceLineKernelMatrix(
        t_pos, s_pos, weight * qp_st_speed_config_.total_time() / t_pos.size());
  }
  return true;
}

// 根据在time时刻各个boundary的类别，确定s的上界或下界
Status QpSplineStGraph::GetSConstraintByTime(const std::vector<const StBoundary*>& boundaries,
                                             const double                          time,
                                             const double                          total_path_s,
                                             double* const                         s_upper_bound,
                                             double* const s_lower_bound) const {
  *s_upper_bound = total_path_s;

  for (const StBoundary* boundary : boundaries) {
    double s_upper = 0.0;
    double s_lower = 0.0;

    if (!boundary->GetUnblockSRange(time, &s_upper, &s_lower)) { continue; }

    // 如果boundary的类型是如下3种，则寻找s的上界，s上界之上可以理解为有障碍物
    // 如follow场景，我们只在意s的上界
    if (boundary->boundary_type() == StBoundary::BoundaryType::STOP ||
        boundary->boundary_type() == StBoundary::BoundaryType::FOLLOW ||
        boundary->boundary_type() == StBoundary::BoundaryType::YIELD) {
      *s_upper_bound = std::fmin(*s_upper_bound, s_upper);
    } else if (boundary->boundary_type() == StBoundary::BoundaryType::OVERTAKE) {
      // 如果boundary的类型是OVERTAKE超车，则寻找s的下界，s下界之下可以理解为有障碍物
      // 即自车要超过障碍物，则自车的s一定要大于障碍物的boundary的upper
      *s_lower_bound = std::fmax(*s_lower_bound, s_lower);
    } else {
      AWARN << "Unhandled boundary type: " << StBoundary::TypeName(boundary->boundary_type());
    }
  }

  return Status::OK();
}

const SpeedData QpSplineStGraph::GetHistorySpeed() const {
  const auto* last_frame = FrameHistory::instance()->Latest();
  if (!last_frame) {
    AWARN << "last frame is empty";
    return SpeedData();
  }
  const ReferenceLineInfo* last_reference_line_info = last_frame->DriveReferenceLineInfo();
  if (!last_reference_line_info) {
    ADEBUG << "last reference line info is empty";
    return SpeedData();
  }
  return last_reference_line_info->speed_data();
}
// 查找针对每个采样时刻t的速度上限
Status QpSplineStGraph::EstimateSpeedUpperBound(const common::TrajectoryPoint& init_point,
                                                const SpeedLimit&              speed_limit,
                                                std::vector<double>* speed_upper_bound) const {
  DCHECK_NOTNULL(speed_upper_bound);

  speed_upper_bound->clear();

  // use v to estimate position: not accurate, but feasible in cyclic
  // processing. We can do the following process multiple times and use
  // previous cycle's results for better estimation.
  const double v               = init_point.v();
  auto         last_speed_data = GetHistorySpeed();

  if (static_cast<double>(t_evaluated_.size() + speed_limit.speed_limit_points().size()) <
      t_evaluated_.size() *
          std::log(static_cast<double>(speed_limit.speed_limit_points().size()))) {
    uint32_t i = 0;
    uint32_t j = 0;

    // 2种不同的方式粗略估计时刻t可能会达到的s，1根据初始速度，2根据上一次规划的速度
    while (i < t_evaluated_.size() && j + 1 < speed_limit.speed_limit_points().size()) {
      double distance = v * t_evaluated_[i];
      if (!last_speed_data.Empty() && distance < last_speed_data.speed_vector().back().s()) {
        SpeedPoint p;
        last_speed_data.EvaluateByTime(t_evaluated_[i], &p);
        distance = p.s();
      }
      // 根据估算的s 求匹配的speed limit
      constexpr double kDistanceEpsilon = 1e-6;
      if (fabs(distance - speed_limit.speed_limit_points()[j].first) < kDistanceEpsilon) {
        speed_upper_bound->push_back(speed_limit.speed_limit_points()[j].second);
        ++i;
      } else if (distance < speed_limit.speed_limit_points()[j].first) {
        ++i;
      } else if (distance <= speed_limit.speed_limit_points()[j + 1].first) {
        speed_upper_bound->push_back(speed_limit.GetSpeedLimitByS(distance));
        ++i;
      } else {
        ++j;
      }
    }

    for (uint32_t k = speed_upper_bound->size(); k < t_evaluated_.size(); ++k) {
      speed_upper_bound->push_back(FLAGS_planning_upper_speed_limit);
      ADEBUG << "speed upper bound:" << speed_upper_bound->back();
    }
  } else {
    auto cmp = [](const std::pair<double, double>& p1, const double s) { return p1.first < s; };

    const auto& speed_limit_points = speed_limit.speed_limit_points();
    for (const double t : t_evaluated_) {
      double s = v * t;
      if (!last_speed_data.Empty() && s < last_speed_data.speed_vector().back().s()) {
        SpeedPoint p;
        last_speed_data.EvaluateByTime(t, &p);
        s = p.s();
      }

      // NOTICE: we are using binary search here based on two assumptions:
      // (1) The s in speed_limit_points increase monotonically.
      // (2) The evaluated_t_.size() << number of speed_limit_points.size()
      //
      // If either of the two assumption is failed, a new algorithm must be
      // used to replace the binary search.

      const auto& it =
          std::lower_bound(speed_limit_points.begin(), speed_limit_points.end(), s, cmp);
      if (it != speed_limit_points.end()) {
        speed_upper_bound->push_back(it->second);
      } else {
        speed_upper_bound->push_back(speed_limit_points.back().second);
      }
    }
  }

  // 如果正在变道，略增大速度上限（0.05）
  if (is_change_lane_) {
    for (uint32_t k = 0; k < t_evaluated_.size(); ++k) {
      speed_upper_bound->at(k) *= (1.0 + FLAGS_change_lane_speed_relax_percentage);
    }
  }

  // clang-format off
  const double kTimeBuffer  = 1.0;
  const double kSpeedBuffer = 0.1;
  // 提高前1s速度上限
  for (uint32_t k = 0;  k < t_evaluated_.size() && t_evaluated_[k] < kTimeBuffer; 
    ++k) {
    speed_upper_bound->at(k) = std::fmax(init_point_.v() 
                             + kSpeedBuffer, speed_upper_bound->at(k));
  }
  // clang-format on
  return Status::OK();
}

}  // namespace planning
}  // namespace apollo
