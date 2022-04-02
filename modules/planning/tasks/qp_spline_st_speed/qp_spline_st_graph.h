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
 * @file qp_spline_st_graph.h
 **/

#ifndef MODULES_PLANNING_TASKS_QP_SPLINE_ST_SPEED_QP_SPLINE_ST_GRAPH_H_
#define MODULES_PLANNING_TASKS_QP_SPLINE_ST_SPEED_QP_SPLINE_ST_GRAPH_H_

#include <memory>
#include <utility>
#include <vector>

#include "modules/common/configs/proto/vehicle_config.pb.h"
#include "modules/common/proto/pnc_point.pb.h"
#include "modules/planning/proto/qp_st_speed_config.pb.h"

#include "modules/common/status/status.h"
#include "modules/common/util/string_util.h"
#include "modules/planning/common/path/path_data.h"
#include "modules/planning/common/planning_util.h"
#include "modules/planning/common/speed/speed_data.h"
#include "modules/planning/common/speed/st_boundary.h"
#include "modules/planning/math/smoothing_spline/spline_1d_generator.h"
#include "modules/planning/tasks/st_graph/st_graph_data.h"

namespace apollo {
namespace planning {

/*
Apollo
该算法的精妙之处就在于，将path和speed分别在SL和ST空间中进行考虑，
使得两者的优化思想非常类似，很巧妙地完成两个维度的求解。但与此同时，
我感觉这也限制了speed优化，对于动态障碍物的处理就不够完备，后面我会单独再写文章详细讲这一块。
*/
class QpSplineStGraph {
 public:
  QpSplineStGraph(Spline1dGenerator*                  spline_generator,
                  const QpStSpeedConfig&              qp_st_speed_config,
                  const apollo::common::VehicleParam& veh_param,
                  const bool                          is_change_lane);

  void SetDebugLogger(planning_internal::STGraphDebug* st_graph_debug);

  common::Status Search(const StGraphData&               st_graph_data,
                        const std::pair<double, double>& accel_bound,
                        const SpeedData&                 reference_speed_data,
                        SpeedData* const                 speed_data);

 private:
  void Init();

  // Add st graph constraint
  common::Status AddConstraint(const common::TrajectoryPoint&        init_point,
                               const SpeedLimit&                     speed_limit,
                               const std::vector<const StBoundary*>& boundaries,
                               const std::pair<double, double>&      accel_bound);

  // Add objective function
  common::Status AddKernel(const std::vector<const StBoundary*>& boundaries,
                           const SpeedLimit&                     speed_limit);
  common::Status Solve();  // solve
  // extract upper lower bound for constraint;
  common::Status GetSConstraintByTime(const std::vector<const StBoundary*>& boundaries,
                                      const double                          time,
                                      const double                          total_path_s,
                                      double* const                         s_upper_bound,
                                      double* const                         s_lower_bound) const;

  // reference line kernel is a constant s line at s = 250m
  common::Status AddCruiseReferenceLineKernel(const double weight);

  // follow line kernel
  common::Status AddFollowReferenceLineKernel(const std::vector<const StBoundary*>& boundaries,
                                              const double                          weight);

  // yield line kernel
  common::Status  AddYieldReferenceLineKernel(const std::vector<const StBoundary*>& boundaries,
                                              const double                          weight);
  const SpeedData GetHistorySpeed() const;
  common::Status  EstimateSpeedUpperBound(const common::TrajectoryPoint& init_point,
                                          const SpeedLimit&              speed_limit,
                                          std::vector<double>*           speed_upper_bound) const;

  bool AddDpStReferenceKernel(const double weight) const;

 private:
  Spline1dGenerator*               spline_generator_ = nullptr;    // solver
  const QpStSpeedConfig            qp_st_speed_config_;            // qp st configuration
  common::TrajectoryPoint          init_point_;                    // initial status
  bool                             is_change_lane_     = false;    // is change lane
  double                           t_knots_resolution_ = 0.0;      // t knots resolution
  std::vector<double>              t_knots_;                       // knots
  double                           t_evaluated_resolution_ = 0.0;  // evaluated t resolution
  std::vector<double>              t_evaluated_;                   // evaluated points
  std::vector<double>              cruise_;                        // reference line kernel
  planning_internal::STGraphDebug* st_graph_debug_ = nullptr;
  // reference st points from dp optimizer
  std::vector<common::SpeedPoint> reference_dp_speed_points_;
};

}  // namespace planning
}  // namespace apollo

#endif  // MODULES_PLANNING_TASKS_QP_SPLINE_ST_SPEED_QP_SPLINE_ST_GRAPH_H_
