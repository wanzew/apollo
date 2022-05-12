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
 * @file trajectory_analyzer.h
 * @brief Defines the TrajectoryAnalyzer class.
 */

#pragma once

#include <vector>

#include "modules/common/proto/pnc_point.pb.h"
#include "modules/planning/proto/planning.pb.h"

#include "modules/common/math/box2d.h"
#include "modules/common/math/vec2d.h"
#include "modules/common/vehicle_state/vehicle_state_provider.h"

/**
 * @namespace apollo::control
 * @brief apollo::control
 */
namespace apollo {
namespace control {

/**
 * @class TrajectoryAnalyzer
 * @brief process point query and conversion related to trajectory
 */
class TrajectoryAnalyzer {
 public:
  /**
   * @brief constructor
   */
  TrajectoryAnalyzer() = default;

  /**
   * @brief constructor
   * @param planning_published_trajectory trajectory data generated by
   * planning module
   */
  TrajectoryAnalyzer(const planning::ADCTrajectory* planning_published_trajectory);

  /**
   * @brief destructor
   */
  ~TrajectoryAnalyzer() = default;

  /**
   * @brief get sequence number of the trajectory
   * @return sequence number.
   */
  unsigned int seq_num() { return seq_num_; }

  /**
   * @brief query a point of trajectory that its absolute time is closest
   * to the give time.
   * @param t absolute time for query
   * @return a point of trajectory
   */
  common::TrajectoryPoint QueryNearestPointByAbsoluteTime(const double t) const;

  /**
   * @brief query a point of trajectory that its relative time is closest
   * to the give time. The time is relative to the first pointof trajectory
   * @param t relative time for query
   * @return a point of trajectory
   */
  common::TrajectoryPoint QueryNearestPointByRelativeTime(const double t) const;

  /**
   * @brief query a point of trajectory that its position is closest
   * to the given position.
   * @param x value of x-coordination in the given position
   * @param y value of y-coordination in the given position
   * @return a point of trajectory
   */
  common::TrajectoryPoint QueryNearestPointByPosition(const double x, const double y) const;

  /**
   * @brief query a point on trajectory that its position is closest
   * to the given position.
   * @param x value of x-coordination in the given position
   * @param y value of y-coordination in the given position
   * @return a point on trajectory, the point may be a point of trajectory
   * or interpolated by two adjacent points of trajectory
   */
  common::PathPoint QueryMatchedPathPoint(const double x, const double y) const;

  /**
   * @brief convert a position with theta and speed to trajectory frame,
   * - longitudinal and lateral direction to the trajectory
   * @param x x-value of the position
   * @param y y-value of the position
   * @param theta heading angle on the position
   * @param v speed on the position
   * @param matched_point matched point on trajectory for the given position
   * @param ptr_s longitudinal distance
   * @param ptr_s_dot longitudinal speed
   * @param ptr_d lateral distance
   * @param ptr_d_dot lateral speed
   */
  void ToTrajectoryFrame(const double             x,
                         const double             y,
                         const double             theta,
                         const double             v,
                         const common::PathPoint& matched_point,
                         double*                  ptr_s,
                         double*                  ptr_s_dot,
                         double*                  ptr_d,
                         double*                  ptr_d_dot) const;

  /**
   * @brief Transform the current trajectory points to the center of mass(COM)
   * of the vehicle, given the distance from rear wheels to the center of mass.
   * @param rear_to_com_distance Distance from rear wheels to
   *        the vehicle's center of mass.
   */
  void TrajectoryTransformToCOM(const double rear_to_com_distance);

  /**
   * @brief Compute the position of center of mass(COM) of the vehicle,
   *        given the distance from rear wheels to the center of mass.
   * @param rear_to_com_distance Distance from rear wheels to
   *        the vehicle's center of mass.
   * @param path_point PathPoint along the published planning trajectory.
   * @return The position of the vehicle's center of mass.
   */
  common::math::Vec2d ComputeCOMPosition(const double             rear_to_com_distance,
                                         const common::PathPoint& path_point) const;

  /**
   * @brief get all points of the trajectory
   * @return a vector of trajectory points
   */
  const std::vector<common::TrajectoryPoint>& trajectory_points() const;

 private:
  common::PathPoint FindMinDistancePoint(const common::TrajectoryPoint& p0,
                                         const common::TrajectoryPoint& p1,
                                         const double                   x,
                                         const double                   y) const;

  std::vector<common::TrajectoryPoint> trajectory_points_;

  double       header_time_ = 0.0;
  unsigned int seq_num_     = 0;
};

}  // namespace control
}  // namespace apollo
