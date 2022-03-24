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
 * @file dp_st_graph.cc
 **/

#include "modules/planning/tasks/dp_st_speed/dp_st_graph.h"

#include <algorithm>
#include <limits>
#include <string>
#include <utility>

#include "modules/common/proto/pnc_point.pb.h"

#include "modules/common/log.h"
#include "modules/common/math/vec2d.h"
#include "modules/planning/common/planning_gflags.h"
#include "modules/planning/common/planning_thread_pool.h"

namespace apollo {
namespace planning {

using apollo::common::ErrorCode;
using apollo::common::SpeedPoint;
using apollo::common::Status;
using apollo::common::VehicleParam;
using apollo::common::math::Vec2d;

usinge StGraphNode = StGraphPoint;

namespace {
constexpr float kInf = std::numeric_limits<float>::infinity();

bool CheckOverlapOnDpStGraph(const std::vector<const StBoundary*>& boundaries,
                             const StGraphNode&                    p1,
                             const StGraphNode&                    p2) {
  const common::math::LineSegment2d seg(p1.point(), p2.point());
  for (const auto* boundary : boundaries) {
    if (boundary->boundary_type() == StBoundary::BoundaryType::KEEP_CLEAR) { continue; }
    if (boundary->HasOverlap(seg)) { return true; }
  }
  return false;
}
}  // namespace

DpStGraph::DpStGraph(const StGraphData&                      st_graph_data,
                     const DpStSpeedConfig&                  dp_config,
                     const std::vector<const PathObstacle*>& obstacles,
                     const common::TrajectoryPoint&          init_point,
                     const SLBoundary&                       adc_sl_boundary)
    : st_graph_data_(st_graph_data)
    , dp_st_speed_config_(dp_config)
    , obstacles_(obstacles)
    , init_point_(init_point)
    , dp_st_cost_(dp_config, obstacles, init_point_)
    , adc_sl_boundary_(adc_sl_boundary) {
  dp_st_speed_config_.set_total_path_length(
      std::fmin(dp_st_speed_config_.total_path_length(), st_graph_data_.path_data_length()));
  unit_s_ =
      dp_st_speed_config_.total_path_length() / (dp_st_speed_config_.matrix_dimension_s() - 1);
  unit_t_ = dp_st_speed_config_.total_time() / (dp_st_speed_config_.matrix_dimension_t() - 1);
}

Status DpStGraph::Search(SpeedData* const speed_data) {
  constexpr float kBounadryEpsilon = 1e-2;
  for (const auto& boundary : st_graph_data_.st_boundaries()) {
    if (boundary->boundary_type() == StBoundary::BoundaryType::KEEP_CLEAR) { continue; }
    if (boundary->IsPointInBoundary({0.0, 0.0}) ||
        (std::fabs(boundary->min_t()) < kBounadryEpsilon &&
         std::fabs(boundary->min_s()) < kBounadryEpsilon)) {
      std::vector<SpeedPoint> speed_profile;
      float                   t = 0.0;
      for (int i = 0; i < dp_st_speed_config_.matrix_dimension_t(); ++i, t += unit_t_) {
        SpeedPoint speed_point;
        speed_point.set_s(0.0);
        speed_point.set_t(t);
        speed_profile.emplace_back(speed_point);
      }
      speed_data->set_speed_vector(speed_profile);
      return Status::OK();
    }
  }

  if (st_graph_data_.st_boundaries().empty()) {
    ADEBUG << "No path obstacles, dp_st_graph output default speed profile.";
    std::vector<SpeedPoint> speed_profile;
    float                   s = 0.0;
    float                   t = 0.0;
    for (int i = 0; i < dp_st_speed_config_.matrix_dimension_t() &&
                    i < dp_st_speed_config_.matrix_dimension_s();
         ++i, t += unit_t_, s += unit_s_) {
      SpeedPoint speed_point;
      speed_point.set_s(s);
      speed_point.set_t(t);
      const float v_default = unit_s_ / unit_t_;
      speed_point.set_v(v_default);
      speed_point.set_a(0.0);
      speed_profile.emplace_back(std::move(speed_point));
    }
    speed_data->set_speed_vector(std::move(speed_profile));
    return Status::OK();
  }

  InitCostTable();
  CalculateTotalCost();
  RetrieveSpeedProfile(speed_data);

  return Status::OK();
}

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
Status DpStGraph::InitCostTable() {
  uint32_t dim_s = dp_st_speed_config_.matrix_dimension_s();  // 150 列
  uint32_t dim_t = dp_st_speed_config_.matrix_dimension_t();  // 8 行
  DCHECK_GT(dim_s, 2);
  DCHECK_GT(dim_t, 2);
  cost_table_  = std::vector<std::vector<StGraphNode>>(  // 8 行 150 列
      dim_t, std::vector<StGraphNode>(dim_s, StGraphNode()));
  float curr_t = 0.0;
  for (uint32_t t = 0; t < cost_table_.size(); ++t, curr_t += unit_t_) {
    auto& cost_table_t = cost_table_[t];
    float curr_s       = 0.0;
    for (uint32_t s = 0; s < cost_table_t.size(); ++s, curr_s += unit_s_) {
      cost_table_t[s].Init(t, s, STPoint(curr_s, curr_t));
    }
  }
  return Status::OK();
}
/*
 *     t  dim_t = 8
 *     |
 *     |
 *   7 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
 *     |
 *   6 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
 *     |
 *   5 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
 *     |
 *   4 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
 *     |
 *   3 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
 *     |
 *   2 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
 *     |
 *   1 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
 *     |
 *   0 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
 *     |_____________________________________________________________________________________
 *     s0   0       20       40        60         80       100        120               149
 *
 *  cost_table_ = std::vector<std::vector<StGraphNode>>(
 *      dim_t, std::vector<StGraphNode>(dim_s, StGraphNode()));
 *
 *
 * s0  &----------@-----------@-----------@-----------@-----------@-----------@-----------@-> s
 * t0 0s          1s          2s          3s          4s          5s          6s          7s
 */

Status DpStGraph::CalculateTotalCost() {
  // 首先定义了两个变量，这两个变量限定的是下一列的s范围。因为我们在初始化CostTable时，
  // 每个时间点上的采样范围都是[0,total_s_]，在实际的计算过程中，并不需要每次都遍历所有的点，
  // 因此，在通过函数CalculateCostAt()计算当前点的cost过后，需要通过函数GetColRange()计算出下一列的s的范围
  // [next_row_lower,next_row_upper]。
  //  CalculateTotalCost()的外层循环是
  // for (size_t t = 0; t < cost_table_.size(); ++t)
  // 内层循环是
  // for (size_t s = next_row_lower; s <= next_row_upper; ++s)
  // 在这两层循环中，调用CalculateCostAt()，在每个时刻t，遍历该时刻所有的采样点，
  // 采样点的s范围为[next_row_lower,next_row_upper]。

  uint32_t next_row_upper = 0;
  uint32_t next_row_lower = 0;
  for (size_t t = 0; t < cost_table_.size(); ++t) {
    // 最高列，即最大加速度情形下所允许的最大距离采样值 default 0
    int upper_row = 0;
    // 最低列，即最大减速度情形下所允许的最小距离采样值 default 149
    int lower_row = cost_table_.back().size() - 1;

    for (uint32_t s = next_row_lower; s <= next_row_upper; ++s) {
      if (FLAGS_enable_multi_thread_in_dp_st_graph) {
        PlanningThreadPool::instance()->Push(std::bind(&DpStGraph::CalculateCostAt, this, t, s));
      } else {
        CalculateCostAt(t, s);
      }
    }
    if (FLAGS_enable_multi_thread_in_dp_st_graph) { PlanningThreadPool::instance()->Synchronize(); }

    /*
     *     t  dim_t = 8
     *     |
     *     |
     *   7 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  #  #  #  #  #  #  #
     *     |
     *   6 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  #  #  #  #  #  #  #  *  *  *  *
     *     |
     *   5 |    *  *  *  *  *  *  *  *  *  *  *  *  #  #  #  #  #  #  #  *  *  *  *  *  *  *  *
     *     |
     *   4 |    *  *  *  *  *  *  *  *  #  #  #  #  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   3 |    *  *  *  *  *  #  #  #  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   2 |    *  *  *  #  #  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   1 |    *  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   0 |    @  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |_____________________________________________________________________________________s
     *     s0   0       20       40        60         80       100        120               149
     *
     */
    // 对于第一、二，直至最后一个时间采样值，循环计算不同距离采样值上的代价
    // 注意：每个时间采样值上的距离采样值数目是不相同的。例如：
    // 第一个时间采样值为起点，在该点上只能有一个距离采样值：0，否则
    // 代价表cost_table_就不正确。正常情况下，第二个时间采样值上的距离采样值数目
    // 会大于1，不然就是匀速前进，玩不下去了
    for (uint32_t s = next_row_lower; s <= next_row_upper; ++s) {
      const auto& cur_node = cost_table_[t][s];
      if (cur_node.total_cost() < std::numeric_limits<float>::infinity()) {
        int h_s = 0;
        int l_s = 0;
        GetColRange(cur_node, &h_s, &l_s);
        upper_row = std::max(upper_row, h_s);
        lower_row = std::min(lower_row, l_s);
      }
    }
    next_row_upper = upper_row;
    next_row_lower = lower_row;
  }

  return Status::OK();
}  // namespace planning

void DpStGraph::GetColRange(const StGraphNode& point, int* next_row_upper, int* next_row_lower) {
  float v0 = 0.0;
  if (!point.pre_point()) {
    v0 = init_point_.v();
  } else {
    v0 = (point.index_s() - point.pre_point()->index_s()) * unit_s_ / unit_t_;
  }

  const int   max_s_size          = cost_table_.back().size() - 1;
  const float speed_coeff         = unit_t_ * unit_t_;
  const float delta_s_upper_bound = v0 * unit_t_ + vehicle_param_.max_acceleration() * speed_coeff;

  *next_row_upper = point.index_s() + static_cast<int>(delta_s_upper_bound / unit_s_);
  if (*next_row_upper >= max_s_size) { *next_row_upper = max_s_size; }

  const float delta_s_lower_bound =
      std::fmax(0.0, v0 * unit_t_ + vehicle_param_.max_deceleration() * speed_coeff);
  *next_row_lower = point.index_s() + static_cast<int>(delta_s_lower_bound / unit_s_);

  if (*next_row_lower > max_s_size) {
    *next_row_lower = max_s_size;
  } else if (*next_row_lower < 0) {
    *next_row_lower = 0;
  }
}

// 通过动态规划的方法进行速度规划就是在CalculateCostAt()中进行的。CalculateCostAt()输入的参数是(t,s)。
// 注意这里的c和r指的是列和行的序号，而不是具体坐标。在Apollo的代码中，表示点的行列号用(t,s)，
// 具体坐标用(curr_t,curr_s)。CalculateCostAt()中的输入参数是[t,s]，表示的是采样点的index。
// 通过遍历cost_table_中所有采样点，在CalculateCostAt()中计算当前点cost_ts的TotalCost。
// 首先，通过GetObstacleCost()计算采样点的obstacle_cost，如果obstacle_cost无穷大，则返回，计算下一个点。
// 通过GetSpatialPotentialCost()计算SpatialPotentialCost，这部分反映的是当前点到终点的距离的cost。
// const auto& cost_init = cost_table_[0][0];
// 定义cost_table_的第一个点为起始点。
// if (t == 0){...}
// if (t == 1){...}
// if (t == 2){...}
// 这里是根据cost_table_的列数（即在采样的时刻的数目）的不同进行分类计算，而不是在判断cost的值。
// 在这里有的博客解析的错误的，请大家注意！！！实际上是分了四类，c==2之后进行的是t>2的计算，只是代码中并没有具体标出来。
// 具体来看：
// if (t == 0) {
//   DCHECK_EQ(s, 0) << "Incorrect. s should be 0 with t = 0. row: " << s;
//   cur_node.SetTotalCost(0.0);//起点的cost设置为0
//   return;
// }
//  当c==0，也就是cost_table_的第一列，t=0时刻，这里只有一个点——起始点cost_table_[0][0]。起始点的TotalCost自然设置为0。
//
void DpStGraph::CalculateCostAt(const uint32_t t, const uint32_t s) {
  auto& cur_node = cost_table_[t][s];
  cur_node.SetObstacleCost(dp_st_cost_.GetObstacleCost(cur_node));

  // 如果障碍物的代价很大，就不再计算别的代价了()
  if (cur_node.obstacle_cost() > std::numeric_limits<float>::max()) { return; }

  const auto& cost_init = cost_table_[0][0];
  // t 等于 0 时， s应该等于0 表示是起点
  if (t == 0) {
    DCHECK_EQ(s, 0) << "Incorrect. s should be 0 with t = 0. row: " << s;
    cur_node.SetTotalCost(0.0);
    return;
    /*
     *     t  dim_t = 8
     *     |
     *     |
     *   7 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  #  #  #  #  #  #  #
     *     |
     *   6 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  #  #  #  #  #  #  #  *  *  *  *
     *     |
     *   5 |    *  *  *  *  *  *  *  *  *  *  *  *  #  #  #  #  #  #  #  *  *  *  *  *  *  *  *
     *     |
     *   4 |    *  *  *  *  *  *  *  *  #  #  #  #  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   3 |    *  *  *  *  *  #  #  #  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   2 |    *  *  *  #  #  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   1 |    *  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   0 |    @  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |_____________________________________________________________________________________s
     *     s0   0       20       40        60         80       100        120               149
     *
     */
  }

  float speed_limit = st_graph_data_.speed_limit().GetSpeedLimitByS(unit_s_ * s);

  if (t == 1) {
    const float acc = (s * unit_s_ / unit_t_ - init_point_.v()) / unit_t_;
    // 如果加速度过大过过小，直接返回，cost默认无穷大
    if (acc < dp_st_speed_config_.max_deceleration() ||
        acc > dp_st_speed_config_.max_acceleration()) {
      return;
    }
    // 如果和起点连线与boundary相交，直接return，cost默认为无穷大
    if (CheckOverlapOnDpStGraph(st_graph_data_.st_boundaries(), cur_node, cost_init)) { return; }

    // 代价叠加
    cur_node.SetTotalCost(cur_node.obstacle_cost() + cost_init.total_cost() +
                          CalculateEdgeCostForSecondCol(s, speed_limit));
    cur_node.SetPrePoint(cost_init);
    return;

    /*
     *     t  dim_t = 8
     *     |
     *     |
     *   7 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  #  #  #  #  #  #  #
     *     |
     *   6 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  #  #  #  #  #  #  #  *  *  *  *
     *     |
     *   5 |    *  *  *  *  *  *  *  *  *  *  *  *  #  #  #  #  #  #  #  *  *  *  *  *  *  *  *
     *     |
     *   4 |    *  *  *  *  *  *  *  *  #  #  #  #  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   3 |    *  *  *  *  *  #  #  #  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   2 |    *  *  *  #  #  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   1 |    *  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   0 |    @  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |_____________________________________________________________________________________s
     *     s0   0       20       40        60         80       100        120               149
     *
     */
  }

  constexpr float kSpeedRangeBuffer = 0.20;
  const uint32_t  max_s_diff        = static_cast<uint32_t>(FLAGS_planning_upper_speed_limit *
                                                    (1 + kSpeedRangeBuffer) * unit_t_ / unit_s_);
  const uint32_t  s_low             = (max_s_diff < s ? s - max_s_diff : 0);
  const auto&     pre_t             = cost_table_[t - 1];

  // https://blog.csdn.net/qq_41324346/article/details/105285029
  // 在完成缩小前一列计算范围的工作之后，计算第三列采样点的TotalCost。遍历前一列[r_low,
  // r]范围内的点，
  // 计算，由起始点开始，经第二列的点，到达当前点的TotalCost，计算方法同上面的方法相同。
  // 不同之处在于，由于当前点可经由第二列不同的点到达，因此，会得到几个TotalCost，从这些TotalCost中找到最小的一个，
  // 并记录下对应的PrePoint，记为当前点的PrePoint。

  if (t == 2) {
    // 寻找前继节点
    for (uint32_t s_pre = s_low; s_pre <= s; ++s_pre) {
      // TODO: 弄懂咋算的
      const float acc = (s * unit_s_ - 2 * s_pre * unit_s_) / (unit_t_ * unit_t_);
      if (acc < dp_st_speed_config_.max_deceleration() ||
          acc > dp_st_speed_config_.max_acceleration()) {
        continue;
      }

      if (CheckOverlapOnDpStGraph(st_graph_data_.st_boundaries(), cur_node, pre_t[s_pre])) {
        continue;
      }

      const float cost = cur_node.obstacle_cost() + pre_t[s_pre].total_cost() +
                         CalculateEdgeCostForThirdCol(s, s_pre, speed_limit);

      if (cost < cur_node.total_cost()) {
        cur_node.SetTotalCost(cost);
        cur_node.SetPrePoint(pre_t[s_pre]);
      }
    }
    return;
    /*
     *     t  dim_t = 8
     *     |
     *     |
     *   7 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  #  #  #  #  #  #  #
     *     |
     *   6 |    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  #  #  #  #  #  #  #  *  *  *  *
     *     |
     *   5 |    *  *  *  *  *  *  *  *  *  *  *  *  #  #  #  #  #  #  #  *  *  *  *  *  *  *  *
     *     |
     *   4 |    *  *  *  *  *  *  *  *  #  #  #  #  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   3 |    *  *  *  *  *  #  #  #  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   2 |    *  *  *  %  #  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   1 |    *  #  #  #  #  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |
     *   0 |    @  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
     *     |_____________________________________________________________________________________s
     *     s0   0       20       40        60         80       100        120               149
     *
     */
  }

  // 当计算进行到第4列及以后，就会执行这部分代码。同样的，在计算之前，已经进行了缩小前一列范围的工作，
  // 这里遍历前一列[r_low, r]范围内的点，计算当前点的cost。这一部分的计算方法与c=2时的是相同的，
  // 之所以再写一遍是计算JerkCost时，这里会多计算一个点，仅此处不同而已。
  // 同样的，在最后会找到当前列TotalCost最小的点，并记录下它的PrePoint。
  // 以上就是CalculateCostAt()的全部内容。大家可以看到，在逐列计算的过程中，
  // 我们已经找到每一列TotalCost最小的点，并且将它与前一列的点关联起来。
  // 只要确定了终点，就可以拔出萝卜带出泥，找到从起点到终点所经过的所有采样点。
  if (t >= 3) {
    for (uint32_t s_pre = s_low; s_pre <= s; ++s_pre) {
      if (std::isinf(pre_t[s_pre].total_cost()) || pre_t[s_pre].pre_point() == nullptr) {
        continue;
      }

      const float curr_a =
          (cur_node.index_s() * unit_s_ + pre_t[s_pre].pre_point()->index_s() * unit_s_ -
           2 * pre_t[s_pre].index_s() * unit_s_) /
          (unit_t_ * unit_t_);
      if (curr_a > vehicle_param_.max_acceleration() ||
          curr_a < vehicle_param_.max_deceleration()) {
        continue;
      }
      if (CheckOverlapOnDpStGraph(st_graph_data_.st_boundaries(), cur_node, pre_t[s_pre])) {
        continue;
      }

      uint32_t           s_prepre           = pre_t[s_pre].pre_point()->index_s();
      const StGraphNode& prepre_graph_point = cost_table_[t - 2][s_prepre];
      if (std::isinf(prepre_graph_point.total_cost())) { continue; }

      if (!prepre_graph_point.pre_point()) { continue; }
      const STPoint& triple_pre_point = prepre_graph_point.pre_point()->point();
      const STPoint& prepre_point     = prepre_graph_point.point();
      const STPoint& pre_point        = pre_t[s_pre].point();
      const STPoint& curr_point       = cur_node.point();
      float          cost =
          cur_node.obstacle_cost() + pre_t[s_pre].total_cost() +
          CalculateEdgeCost(triple_pre_point, prepre_point, pre_point, curr_point, speed_limit);

      if (cost < cur_node.total_cost()) {
        cur_node.SetTotalCost(cost);
        cur_node.SetPrePoint(pre_t[s_pre]);
      }
    }
  }
}

Status DpStGraph::RetrieveSpeedProfile(SpeedData* const speed_data) {
  float              min_cost       = std::numeric_limits<float>::infinity();
  const StGraphNode* best_end_point = nullptr;
  /*
   *     t  dim_t = 8
   *     |
   *   7 |    ****************************************************[************************@**]
   *     |
   *   6 |    *******************************************[*********************@*****]*********
   *     |
   *   5 |    ******************************[************************@******]******************
   *     |
   *   4 |    *******************[*********************@*****]*********************************
   *     |
   *   3 |    ***********[*******************@****]********************************************
   *     |
   *   2 |    *******[**********@**]***********************************************************
   *     |
   *   1 |    **[*******@***]******************************************************************
   *     |
   *   0 |    @********************************************************************************
   *     |_____________________________________________________________________________________
   *     s0   0       20       40        60         80        100         120               149
   *
   */
  // 找到t=7时cost最小的节点
  // 经过CalculateCostAt()的计算，我们得到了最后一列所有点的TotalCost，通过这个for循环，
  // 我们可以找到最后一列TotalCost最小的采样点，将这个点暂定为最佳的终点，
  // 记为best_end_point，并记下它的TotalCost为min_cost。
  for (const StGraphNode& cur_node : cost_table_.back()) {
    if (!std::isinf(cur_node.total_cost()) && cur_node.total_cost() < min_cost) {
      best_end_point = &cur_node;
      min_cost       = cur_node.total_cost();
    }
  }

  // 对于cost_table_中的每一行，即第一个、第二个、...、最后一个时间采样值上的代价值数组，
  // 其最后一个元素存储的是本级时间采样值上的最小代价节点。
  // 将这些节点与现有最优终点 best_end_point 比较，
  // 不断更新最小代价值min_cost和最优终点best_end_point，
  // 直至找到全局最优终点
  for (const auto& row : cost_table_) {
    const StGraphNode& cur_node = row.back();
    if (!std::isinf(cur_node.total_cost()) && cur_node.total_cost() < min_cost) {
      best_end_point = &cur_node;
      min_cost       = cur_node.total_cost();
    }
  }

  if (best_end_point == nullptr) {
    const std::string msg = "Fail to find the best feasible trajectory.";
    AERROR << msg;
    return Status(ErrorCode::PLANNING_ERROR, msg);
  }

  // 设置最优终点的速度数据，并顺藤摸瓜找出其连接的倒数第二个、倒数第三个直到第一个时间节点
  // 分别设置这些时间节点的速度数据
  // 这里需要插一句，为什么不直接将最后一列TotalCost最小的点定位最佳终点呢？
  // 因为采样时间是我们估计的一个大概的时间，在这个时间之前，动态规划可能已经规划到了路径的终点，
  // 只是因为还没有计算到cost_table_的最后一列，才一直进行计算。因此，我们在判断为最后一列之后，
  // 需要继续判断每一列的最后一个点，是不是已经到达了路径的终点。
  std::vector<SpeedPoint> speed_profile;
  const StGraphNode*      cur_node = best_end_point;
  while (cur_node != nullptr) {
    SpeedPoint speed_point;
    speed_point.set_s(cur_node->point().s());
    speed_point.set_t(cur_node->point().t());
    speed_profile.emplace_back(speed_point);
    cur_node = cur_node->pre_point();
  }
  std::reverse(speed_profile.begin(), speed_profile.end());

  constexpr float kEpsilon = std::numeric_limits<float>::epsilon();
  if (speed_profile.front().t() > kEpsilon || speed_profile.front().s() > kEpsilon) {
    const std::string msg = "Fail to retrieve speed profile.";
    AERROR << msg;
    return Status(ErrorCode::PLANNING_ERROR, msg);
  }
  speed_data->set_speed_vector(speed_profile);
  return Status::OK();
}

float DpStGraph::CalculateEdgeCost(const STPoint& first,
                                   const STPoint& second,
                                   const STPoint& third,
                                   const STPoint& forth,
                                   const float    speed_limit) {
  return dp_st_cost_.GetSpeedCost(third, forth, speed_limit) +
         dp_st_cost_.GetAccelCostByThreePoints(second, third, forth) +
         dp_st_cost_.GetJerkCostByFourPoints(first, second, third, forth);
}

float DpStGraph::CalculateEdgeCostForSecondCol(const uint32_t row, const float speed_limit) {
  float          init_speed = init_point_.v();
  float          init_acc   = init_point_.a();
  const STPoint& pre_point  = cost_table_[0][0].point();
  const STPoint& curr_point = cost_table_[1][row].point();
  return dp_st_cost_.GetSpeedCost(pre_point, curr_point, speed_limit) +
         dp_st_cost_.GetAccelCostByTwoPoints(init_speed, pre_point, curr_point) +
         dp_st_cost_.GetJerkCostByTwoPoints(init_speed, init_acc, pre_point, curr_point);
}

float DpStGraph::CalculateEdgeCostForThirdCol(const uint32_t curr_row,
                                              const uint32_t pre_row,
                                              const float    speed_limit) {
  float          init_speed = init_point_.v();
  const STPoint& first      = cost_table_[0][0].point();
  const STPoint& second     = cost_table_[1][pre_row].point();
  const STPoint& third      = cost_table_[2][curr_row].point();
  return dp_st_cost_.GetSpeedCost(second, third, speed_limit) +
         dp_st_cost_.GetAccelCostByThreePoints(first, second, third) +
         dp_st_cost_.GetJerkCostByThreePoints(init_speed, first, second, third);
}

}  // namespace planning
}  // namespace apollo
