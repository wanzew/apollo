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
 * @file qp_spline_path_generator.cc
 **/
#include "modules/planning/tasks/qp_spline_path/qp_spline_path_generator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "modules/common/proto/pnc_point.pb.h"

#include "modules/common/log.h"
#include "modules/common/macro.h"
#include "modules/common/math/cartesian_frenet_conversion.h"
#include "modules/common/util/string_util.h"
#include "modules/common/util/util.h"
#include "modules/planning/common/planning_gflags.h"
#include "modules/planning/common/planning_util.h"

namespace apollo {
namespace planning {

namespace {
double GetLaneChangeLateralShift(const double v) {
  const double l0    = 2.0;    // shift at v = 0 m/s
  const double v_ref = 20.11;  // reference speed: 45mph = 20.11 m/s
  const double l_ref = 1.4;
  const double b     = l0;
  const double a     = (l_ref - b) / v_ref;
  return a * v + b;
}
}  // namespace

using apollo::common::math::CartesianFrenetConverter;
using apollo::common::math::Vec2d;

QpSplinePathGenerator::QpSplinePathGenerator(Spline1dGenerator*        spline_generator,
                                             const ReferenceLine&      reference_line,
                                             const QpSplinePathConfig& qp_spline_path_config,
                                             const SLBoundary&         adc_sl_boundary)
    : spline_generator_(spline_generator)
    , reference_line_(reference_line)
    , qp_spline_path_config_(qp_spline_path_config)
    , adc_sl_boundary_(adc_sl_boundary) {
  CHECK_GE(qp_spline_path_config_.regularization_weight(), 0.0)
      << "regularization_weight should NOT be negative.";
  CHECK_GE(qp_spline_path_config_.derivative_weight(), 0.0)
      << "derivative_weight should NOT be negative.";
  CHECK_GE(qp_spline_path_config_.second_derivative_weight(), 0.0)
      << "second_derivative_weight should NOT be negative.";
  CHECK_GE(qp_spline_path_config_.third_derivative_weight(), 0.0)
      << "third_derivative_weight should NOT be negative.";
}

void QpSplinePathGenerator::SetDebugLogger(apollo::planning_internal::Debug* debug) {
  planning_debug_ = debug;
}

void QpSplinePathGenerator::SetChangeLane(bool is_change_lane_path) {
  is_change_lane_path_ = is_change_lane_path;
}

// 主要有以下关键步骤：
// 初始化样条曲线 InitSpline()，
// 初始化用来计算自车轨迹横向约束的类 QpFrenetFrame::Init()，
// 添加约束 AddConstraint()，
// 添加目标函数 AddKernel()，
// 优化问题求解器 Solve()，
// 最后将求得的轨迹点封装成 DiscretizedPath。
bool QpSplinePathGenerator::Generate(const std::vector<const PathObstacle*>& path_obstacles,
                                     const SpeedData&                        speed_data,
                                     const common::TrajectoryPoint&          init_point,
                                     const double                            boundary_extension,
                                     bool                                    is_final_attempt,
                                     PathData* const                         path_data) {
  ADEBUG << "Init point: " << init_point.DebugString();
  init_trajectory_point_ = init_point;

  const auto& path_data_history = path_data->path_data_history();
  if (!path_data_history.empty()) { last_discretized_path_ = &path_data_history.back().first; }

  // calc init point
  CalculateFrenetPoint(init_point, &init_frenet_point_);

  if (is_change_lane_path_) { ref_l_ = init_frenet_point_.l(); }
  double start_s = init_frenet_point_.s();

  constexpr double kDefaultPathLength = 50.0;
  double           end_s              = std::fmin(
      init_frenet_point_.s() +
          std::fmax(kDefaultPathLength, init_trajectory_point_.v() * FLAGS_look_forward_time_sec),
      reference_line_.Length());

  constexpr double kMinPathLength = 1.0e-6;
  if (start_s + kMinPathLength > end_s) {
    AERROR << "Path length is too small. Path start_s: " << start_s << ", end_s: " << end_s;
    return false;
  }

  // InitSpline()主要完成对spline segment和纵向s采样点的初始化。
  InitSpline(start_s, end_s);

  QpFrenetFrame qp_frenet_frame(reference_line_, speed_data, init_frenet_point_,
                                qp_spline_path_config_.time_resolution(), evaluated_s_);

  // QpFrenetFrame::Init() 根据已知条件，计算对自车横向轨迹的约束。
  // 1. HDMap和道路参考线对自车横向轨迹的约束
  // 2. 静态动态障碍物及相应的decision对自车横向轨迹的约束
  // initialize qp frenet frame
  qp_frenet_frame.Init(path_obstacles);
  // debug signal
  qp_frenet_frame.LogQpBound(planning_debug_);

  // setup pss path constraint
  AddConstraint(qp_frenet_frame, boundary_extension);

  AddKernel();

  bool is_solved = Solve();

  if (!is_solved && !is_final_attempt) { return false; }

  if (!is_solved) {
    AERROR << "Fail to solve qp_spline_path. Use reference line as qp_path "
              "output.";
    util::DumpPlanningContext();
  }
  ADEBUG << common::util::StrCat("Spline dl:", init_frenet_point_.dl(),
                                 ", ddl:", init_frenet_point_.ddl());

  // extract data
  const Spline1d&                spline = spline_generator_->spline();
  std::vector<common::PathPoint> path_points;

  ReferencePoint ref_point = reference_line_.GetReferencePoint(start_s);
  Vec2d          xy_point  = CartesianFrenetConverter::CalculateCartesianPoint(
      ref_point.heading(), Vec2d(ref_point.x(), ref_point.y()), init_frenet_point_.l());

  const auto xy_diff = xy_point - Vec2d(init_point.path_point().x(), init_point.path_point().y());

  double           s_resolution = (end_s - start_s) / qp_spline_path_config_.num_output();
  constexpr double kEpsilon     = 1e-6;
  for (double s = start_s; s + kEpsilon < end_s; s += s_resolution) {
    double l = init_frenet_point_.l();
    if (is_solved) { l = spline(s) + ref_l_; }
    if (planning_debug_ && planning_debug_->planning_data().sl_frame().size() >= 1) {
      auto sl_point =
          planning_debug_->mutable_planning_data()->mutable_sl_frame(0)->mutable_sl_path()->Add();
      sl_point->set_l(l);
      sl_point->set_s(s);
    }
    double dl  = 0.0;
    double ddl = 0.0;
    if (is_solved) {
      dl  = spline.Derivative(s);
      ddl = spline.SecondOrderDerivative(s);
    }
    ReferencePoint ref_point     = reference_line_.GetReferencePoint(s);
    Vec2d          curr_xy_point = CartesianFrenetConverter::CalculateCartesianPoint(
                              ref_point.heading(), Vec2d(ref_point.x(), ref_point.y()), l) -
                          xy_diff;
    double theta =
        CartesianFrenetConverter::CalculateTheta(ref_point.heading(), ref_point.kappa(), l, dl);
    double kappa =
        CartesianFrenetConverter::CalculateKappa(ref_point.kappa(), ref_point.dkappa(), l, dl, ddl);

    common::PathPoint path_point = common::util::MakePathPoint(curr_xy_point.x(), curr_xy_point.y(),
                                                               0.0, theta, kappa, 0.0, 0.0);
    if (!path_points.empty()) {
      double distance = common::util::DistanceXY(path_points.back(), path_point);
      path_point.set_s(path_points.back().s() + distance);
      if (distance > 1e-4) {
        path_point.set_dkappa((kappa - path_points.back().kappa()) / distance);
      }
    }

    if (path_point.s() > end_s) { break; }
    path_points.emplace_back(std::move(path_point));
  }
  path_data->SetReferenceLine(&reference_line_);
  path_data->SetDiscretizedPath(DiscretizedPath(path_points));
  return true;
}

bool QpSplinePathGenerator::CalculateFrenetPoint(
    const common::TrajectoryPoint& traj_point, common::FrenetFramePoint* const frenet_frame_point) {
  common::SLPoint sl_point;
  if (!reference_line_.XYToSL({traj_point.path_point().x(), traj_point.path_point().y()},
                              &sl_point)) {
    return false;
  }
  frenet_frame_point->set_s(sl_point.s());
  frenet_frame_point->set_l(sl_point.l());

  const double theta = traj_point.path_point().theta();
  const double kappa = traj_point.path_point().kappa();
  const double l     = frenet_frame_point->l();

  ReferencePoint ref_point;
  ref_point = reference_line_.GetReferencePoint(frenet_frame_point->s());

  const double theta_ref  = ref_point.heading();
  const double kappa_ref  = ref_point.kappa();
  const double dkappa_ref = ref_point.dkappa();

  const double dl =
      CartesianFrenetConverter::CalculateLateralDerivative(theta_ref, theta, l, kappa_ref);
  const double ddl = CartesianFrenetConverter::CalculateSecondOrderLateralDerivative(
      theta_ref, theta, kappa_ref, kappa, dkappa_ref, l);
  frenet_frame_point->set_dl(dl);
  frenet_frame_point->set_ddl(ddl);
  return true;
}

// 将纵向区间[start_s, end_s] 按照qp_spline_path_config_.max_spline_length()和
// qp_spline_path_config_.max_constraint_interval()的设置，均匀分割后存入knots_和evaluated_s_，
// 每一段都会对应一条qp_spline_path_config_.spline_order()次多项式曲线。
// 在这里，我认为knots_和evaluated_s_
// 这2个记录s轴采样点的vector是应该完全相同的。不明白代码中为何有2种定义？
bool QpSplinePathGenerator::InitSpline(const double start_s, const double end_s) {
  uint32_t number_of_spline =
      static_cast<uint32_t>((end_s - start_s) / qp_spline_path_config_.max_spline_length() + 1.0);
  number_of_spline = std::max(1u, number_of_spline);
  common::util::uniform_slice(start_s, end_s, number_of_spline, &knots_);

  // spawn a new spline generator
  // 产生number_of_spline条order阶spline
  spline_generator_->Reset(knots_, qp_spline_path_config_.spline_order());

  // set evaluated_s_
  uint32_t constraint_num =
      (end_s - start_s) / qp_spline_path_config_.max_constraint_interval() + 1;
  common::util::uniform_slice(start_s, end_s, constraint_num - 1, &evaluated_s_);
  return (knots_.size() > 1) && !evaluated_s_.empty();
  //难道 max_spline_length 和 max_constraint_interval 可以不相等吗？
  //难道 constraint_num 和 number_of_spline 可以不相等吗？
}

bool QpSplinePathGenerator::AddConstraint(const QpFrenetFrame& qp_frenet_frame,
                                          const double         boundary_extension) {
  Spline1dConstraint* spline_constraint = spline_generator_->mutable_spline_constraint();

  // curve个数 * curve多项式参数个数
  const int        dim         = (knots_.size() - 1) * (qp_spline_path_config_.spline_order() + 1);
  constexpr double param_range = 1e-4;

  //循环（curve个数）次，这是添加什么约束？待后面看看发挥什么作用？
  for (int i = qp_spline_path_config_.spline_order(); i < dim;
       i += qp_spline_path_config_.spline_order() + 1) {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(1, dim);
    Eigen::MatrixXd bd  = Eigen::MatrixXd::Zero(1, 1);
    mat(0, i)           = -1;
    bd(0, 0)            = -param_range;
    spline_constraint->AddInequalityConstraint(mat, bd);
    mat(0, i) = 1;
    bd(0, 0)  = -param_range;
    spline_constraint->AddInequalityConstraint(mat, bd);
  }

  // add init status constraint, equality constraint
  //这里3个函数其实是添加了0，1，2阶导不等式约束 sample points for boundary
  const double kBoundaryEpsilon = 1e-4;
  spline_constraint->AddPointConstraintInRange(init_frenet_point_.s(),
                                               init_frenet_point_.l() - ref_l_, kBoundaryEpsilon);
  spline_constraint->AddPointDerivativeConstraintInRange(init_frenet_point_.s(),
                                                         init_frenet_point_.dl(), kBoundaryEpsilon);

  // path二阶导其实是曲线的曲率kappa，U形急转弯不对kappa添加约束
  if (init_trajectory_point_.v() > qp_spline_path_config_.uturn_speed_limit()) {
    spline_constraint->AddPointSecondDerivativeConstraintInRange(
        init_frenet_point_.s(), init_frenet_point_.ddl(), kBoundaryEpsilon);
  }

  // add end point constraint, equality constraint
  // lat_shift看不懂，为什么要与ref_l_反向？
  //推测：换道时有参考线切换，自车横向坐标在目标车道参考线坐标系下是ref_l_，
  //目标点应该在目标车道参考线上采样，所以目标点在当前车道参考线坐标系下
  //横向坐标大概是fabs(ref_l_)，而正负一定相反
  //如果不变道，ref_l_=0，如果变道，ref_l_=init_frenet_point_.l();
  double lat_shift = -ref_l_;
  if (is_change_lane_path_) {
    double lane_change_lateral_shift = GetLaneChangeLateralShift(init_trajectory_point_.v());
    lat_shift = std::copysign(std::fmin(std::fabs(ref_l_), lane_change_lateral_shift), -ref_l_);
  }

  const double     kEndPointBoundaryEpsilon = 1e-2;
  constexpr double kReservedDistance        = 20.0;
  const double     target_s = std::fmin(qp_spline_path_config_.point_constraint_s_position(),
                                    kReservedDistance + init_frenet_point_.s() +
                                        init_trajectory_point_.v() * FLAGS_look_forward_time_sec);
  spline_constraint->AddPointConstraintInRange(target_s, lat_shift, kEndPointBoundaryEpsilon);
  //如果是变道，则不约束一阶导即方向角，如果不是变道，则目标点的方向=0，即沿s轴
  if (!is_change_lane_path_) {
    spline_constraint->AddPointDerivativeConstraintInRange(evaluated_s_.back(), 0.0,
                                                           kEndPointBoundaryEpsilon);
  }

  // add first derivative bound to improve lane change smoothness
  // dl_bound是0.1，arctan(0.1)=5.71度，这角度会不会太小了？
  std::vector<double> dl_lower_bound(evaluated_s_.size(), -FLAGS_dl_bound);
  std::vector<double> dl_upper_bound(evaluated_s_.size(), FLAGS_dl_bound);

  // 添加不等式约束
  // add derivative boundary
  spline_constraint->AddDerivativeBoundary(evaluated_s_, dl_lower_bound, dl_upper_bound);

  // kappa bound is based on the inequality:
  // kappa = d(phi)/ds <= d(phi)/dx = d2y/dx2
  // kappa = d(phi)/ds <= d(phi)/dx = d2y/dx2，此式如何得来？
  // add second derivative boundary
  std::vector<double> kappa_lower_bound(evaluated_s_.size(), -FLAGS_kappa_bound);
  std::vector<double> kappa_upper_bound(evaluated_s_.size(), FLAGS_kappa_bound);
  spline_constraint->AddSecondDerivativeBoundary(evaluated_s_, kappa_lower_bound,
                                                 kappa_upper_bound);

  // dkappa = d(kappa) / ds <= d3y/dx3
  std::vector<double> dkappa_lower_bound(evaluated_s_.size(), -FLAGS_dkappa_bound);
  std::vector<double> dkappa_upper_bound(evaluated_s_.size(), FLAGS_dkappa_bound);

  // add third derivative boundary
  spline_constraint->AddThirdDerivativeBoundary(evaluated_s_, dkappa_lower_bound,
                                                dkappa_upper_bound);

  // add map bound constraint
  double lateral_buf = boundary_extension;
  if (is_change_lane_path_) { lateral_buf = qp_spline_path_config_.cross_lane_lateral_extension(); }
  std::vector<double> boundary_low;
  std::vector<double> boundary_high;

  for (uint32_t i = 0; i < evaluated_s_.size(); ++i) {
    auto road_boundary        = qp_frenet_frame.GetMapBound().at(i);
    auto static_obs_boundary  = qp_frenet_frame.GetStaticObstacleBound().at(i);
    auto dynamic_obs_boundary = qp_frenet_frame.GetDynamicObstacleBound().at(i);

    if (evaluated_s_.at(i) - evaluated_s_.front() <
        qp_spline_path_config_.cross_lane_longitudinal_extension()) {
      //扩宽边界约束
      //右边界
      road_boundary.first = std::fmin(road_boundary.first, init_frenet_point_.l() - lateral_buf);
      //左边界
      road_boundary.second = std::fmax(road_boundary.second, init_frenet_point_.l() + lateral_buf);
    }
    boundary_low.emplace_back(common::util::MaxElement(std::vector<double>{
        road_boundary.first, static_obs_boundary.first, dynamic_obs_boundary.first}));
    boundary_high.emplace_back(common::util::MinElement(std::vector<double>{
        road_boundary.second, static_obs_boundary.second, dynamic_obs_boundary.second}));
  }

  //不等式约束
  const double start_l = ref_l_;

  //求横向相对坐标
  std::for_each(boundary_low.begin(), boundary_low.end(), [start_l](double& d) { d -= start_l; });
  std::for_each(boundary_high.begin(), boundary_high.end(), [start_l](double& d) { d -= start_l; });

  // Add boundary constraint
  spline_constraint->AddBoundary(evaluated_s_, boundary_low, boundary_high);

  // add spline joint third derivative constraint 等式约束
  spline_constraint->AddThirdDerivativeSmoothConstraint();

  return true;
}

void QpSplinePathGenerator::AddHistoryPathKernel() {
  if (last_discretized_path_ == nullptr) { return; }

  PathData last_path_data;
  last_path_data.SetReferenceLine(&reference_line_);
  last_path_data.SetDiscretizedPath(*last_discretized_path_);

  std::vector<double> history_s;
  std::vector<double> histroy_l;

  for (uint32_t i = 0; i < last_path_data.frenet_frame_path().NumOfPoints(); ++i) {
    const auto p = last_path_data.frenet_frame_path().PointAt(i);
    history_s.push_back(p.s());
    histroy_l.push_back(p.l() - ref_l_);
  }

  Spline1dKernel* spline_kernel = spline_generator_->mutable_spline_kernel();
  spline_kernel->AddReferenceLineKernelMatrix(history_s, histroy_l,
                                              qp_spline_path_config_.history_path_weight());
}

// AddKernel()中，如果处在变道的过程中，则添加一项额外的cost，并且只作用于第一段spline。
// AddDerivativeKernelMatrixForSplineK()与AddDerivativeKernelMatrix()如出一辙，
// 唯一的区别是因为针对特定段的spline，没有对knots的循环。
// 最底层同样都是调用SplineSegKernel::Instance()->NthDerivativeKernel()实现的。
void QpSplinePathGenerator::AddKernel() {
  Spline1dKernel* spline_kernel = spline_generator_->mutable_spline_kernel();

  //只有在!is_change_lane_path_时才会进入，而!is_change_lane_path_会导致ref_l_=0
  if (init_trajectory_point_.v() < qp_spline_path_config_.uturn_speed_limit() &&
      !is_change_lane_path_ && qp_spline_path_config_.reference_line_weight() > 0.0) {
    std::vector<double> ref_l(evaluated_s_.size(), -ref_l_);

    // 添加和参考线相关的kernel，推测应该是使轨迹尽可能贴近参考线，但我在代码中
    // 没有理解这一点。而且每次ref_l元素全是0，会导致函数内offset_全是0，意义何在？
    // 参考线或下面的历史轨迹，可以是DP输出的粗糙规划结果
    spline_kernel->AddReferenceLineKernelMatrix(evaluated_s_, ref_l,
                                                qp_spline_path_config_.reference_line_weight());
  }

  //添加和历史轨迹相关的kernel，推测应该是使轨迹尽可能近似延续之前的轨迹，
  //但我在代码中没有理解这一点
  if (qp_spline_path_config_.history_path_weight() > 0.0) { AddHistoryPathKernel(); }

  if (qp_spline_path_config_.regularization_weight() > 0.0) {
    // 添加正则项
    spline_kernel->AddRegularization(qp_spline_path_config_.regularization_weight());
  }

  if (qp_spline_path_config_.derivative_weight() > 0.0) {
    // 添加1阶导kernel，(f'(x))^2积分
    spline_kernel->AddDerivativeKernelMatrix(qp_spline_path_config_.derivative_weight());
    if (std::fabs(init_frenet_point_.l()) > qp_spline_path_config_.lane_change_mid_l()) {
      spline_kernel->AddDerivativeKernelMatrixForSplineK(
          0, qp_spline_path_config_.first_spline_weight_factor() *
                 qp_spline_path_config_.derivative_weight());
    }
  }

  if (qp_spline_path_config_.second_derivative_weight() > 0.0) {
    // 添加2阶导kernel，(f''(x))^2积分
    spline_kernel->AddSecondOrderDerivativeMatrix(
        qp_spline_path_config_.second_derivative_weight());
    if (std::fabs(init_frenet_point_.l()) > qp_spline_path_config_.lane_change_mid_l()) {
      spline_kernel->AddSecondOrderDerivativeMatrixForSplineK(
          0, qp_spline_path_config_.first_spline_weight_factor() *
                 qp_spline_path_config_.second_derivative_weight());
    }
  }

  if (qp_spline_path_config_.third_derivative_weight() > 0.0) {
    // 添加3阶导kernel，(f'''(x))^2积分
    spline_kernel->AddThirdOrderDerivativeMatrix(qp_spline_path_config_.third_derivative_weight());
    if (std::fabs(init_frenet_point_.l()) > qp_spline_path_config_.lane_change_mid_l()) {
      spline_kernel->AddThirdOrderDerivativeMatrixForSplineK(
          0, qp_spline_path_config_.first_spline_weight_factor() *
                 qp_spline_path_config_.third_derivative_weight());
    }
  }
}

bool QpSplinePathGenerator::Solve() {
  // solve the qp problem in spline path generator
  spline_generator_->Solve();
  return true;
}

}  // namespace planning
}  // namespace apollo
