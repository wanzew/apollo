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

#include "modules/planning/common/speed/st_boundary.h"

#include <algorithm>
#include <cmath>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "modules/common/log.h"

namespace apollo {
namespace planning {

using apollo::common::math::Box2d;
using apollo::common::math::Vec2d;

TEST(StBoundaryTest, basic_test) {
  /* s
   * ^
   * |
   * |
   * |
   * |
   * |
   * upper_point[0]      upper_point[1]
   * *--------------------------*
   * |                          |
   * |                          |
   * |                          |
   * *--------------------------*------------------------------------> t
   * lower_point[0]      lower_point[1]
   */

  std::vector<STPoint> upper_points;
  std::vector<STPoint> lower_points;
  // STPoint(const double s, const double t)
  lower_points.emplace_back(0.0, 0.0);
  lower_points.emplace_back(0.0, 10.0);
  upper_points.emplace_back(5.0, 0.0);
  upper_points.emplace_back(5.0, 10.0);

  point_pairs.emplace_back(lower_points[0], upper_points[0]);
  point_pairs.emplace_back(lower_points[1], upper_points[1]);

  std::vector<std::pair<STPoint, STPoint>> point_pairs;
  // 每个pair包含两个点(t相同)，point_pairs 包含两个或多个元素元素
  // 第一个元素包含的是 boundary(菱形)的左边的两个点(小t)
  // 第二个元素包含的是 boundary(菱形)的右边的两个点(大t)
  point_pairs.emplace_back(STPoint(1.0, 0.0), STPoint(5.0, 0.0));
  point_pairs.emplace_back(STPoint(1.0, 10.0), STPoint(5.0, 10.0));

  StBoundary boundary(point_pairs);

  EXPECT_EQ(boundary.id(), "");
  EXPECT_EQ(boundary.boundary_type(), StBoundary::BoundaryType::UNKNOWN);
  EXPECT_FLOAT_EQ(0.0, boundary.min_s());
  EXPECT_FLOAT_EQ(5.0, boundary.max_s());
  EXPECT_FLOAT_EQ(0.0, boundary.min_t());
  EXPECT_FLOAT_EQ(10.0, boundary.max_t());
}

TEST(StBoundaryTest, boundary_range) {
  std::vector<std::pair<STPoint, STPoint>> point_pairs;

  /* s
   * ^
   * |
   * |(5.0, 0)               (5.0, 10.0)
   * *--------------------------*
   * |                          |
   * |                          |
   * |                          |
   * |                          |
   * |                          |
   * *--------------------------*
   * |(1.0, 0)               (1.0, 10.0)
   * ---------------------------------------------------------------> t
   *
   */
  // 每一个pair的两个点，对应的时间t相同
  point_pairs.emplace_back(STPoint(1.0, 0.0), STPoint(5.0, 0.0));
  point_pairs.emplace_back(STPoint(1.0, 10.0), STPoint(5.0, 10.0));

  StBoundary boundary(point_pairs);

  boundary.SetBoundaryType(StBoundary::BoundaryType::YIELD);
  double       t  = -10.0;
  const double dt = 0.01;
  while (t < 10.0) {
    double low  = 0.0;
    double high = 0.0;
    if (t < 0.0) {
      EXPECT_TRUE(boundary.GetUnblockSRange(t, &high, &low));
      EXPECT_DOUBLE_EQ(low, 0.0);
      EXPECT_DOUBLE_EQ(high, 200.0);
      EXPECT_FALSE(boundary.GetBoundarySRange(t, &high, &low));
    } else {
      EXPECT_TRUE(boundary.GetUnblockSRange(t, &high, &low));
      EXPECT_DOUBLE_EQ(low, 0.0);
      EXPECT_DOUBLE_EQ(high, 1.0);

      EXPECT_TRUE(boundary.GetBoundarySRange(t, &high, &low));
      EXPECT_DOUBLE_EQ(low, 1.0);
      EXPECT_DOUBLE_EQ(high, 5.0);
    }
    t += dt;
  }
}

TEST(StBoundaryTest, get_index_range) {
  std::vector<STPoint> upper_points;
  std::vector<STPoint> lower_points;

  /* s
   * ^
   * |
   * |(5.0, 0)               (5.0, 10.0)
   * *--------------------------*
   * |                          |
   * |                          |
   * |                          |
   * |                          |
   * |                          |
   * *--------------------------*
   * |(1.0, 0)               (1.0, 10.0)
   * ---------------------------------------------------------------> t
   *
   */
  std::vector<std::pair<STPoint, STPoint>> point_pairs;

  lower_points.emplace_back(43.000164837720789, -517957.08587679861);
  lower_points.emplace_back(46.100164825451913, -517955.58587660792);

  upper_points.emplace_back(52.200164801309178, -517957.08587679861);
  upper_points.emplace_back(55.6001647283625, -517955.58587660792);

  StBoundary boundary(point_pairs);

  size_t left_index  = 0;
  size_t right_index = 0;
  /* s
   * ^
   * |                          upper_points[1]
   * |                              *
   * |
   * |
   * |
   * |   upper_points[0]            *
   * |        *                 lower_points[1]
   * |
   * |                          -517955.58587660792
   * |
   * |        *
   * |   lower_points[0]
   * |  -517957.08587679861
   * |
   * |
   *  --------------------------------------------------------------> t
   *
   */

  EXPECT_TRUE(boundary.GetIndexRange(lower_points, -517957.08587679861, &left_index, &right_index));
  EXPECT_EQ(left_index, 0);
  EXPECT_EQ(right_index, 0);

  EXPECT_TRUE(boundary.GetIndexRange(lower_points, -517955.58587660792, &left_index, &right_index));
  EXPECT_EQ(left_index, 0);
  EXPECT_EQ(right_index, 1);

  EXPECT_TRUE(boundary.GetIndexRange(lower_points, -517955.58587660792, &left_index, &right_index));
  EXPECT_EQ(left_index, 0);
  EXPECT_EQ(right_index, 1);

  EXPECT_FALSE(boundary.GetIndexRange(lower_points, 0.0, &left_index, &right_index));
}

TEST(StBoundaryTest, remove_redundant_points) {
  std::vector<std::pair<STPoint, STPoint>> point_pairs;
  point_pairs.emplace_back(STPoint(0.0, 0.0), STPoint(1.0, 0.0));
  point_pairs.emplace_back(STPoint(0.1, 0.2), STPoint(1.1, 0.2));
  point_pairs.emplace_back(STPoint(0.2, 0.3), STPoint(1.2, 0.3));
  point_pairs.emplace_back(STPoint(0.3, 0.4), STPoint(1.3, 0.4));
  point_pairs.emplace_back(STPoint(0.4, 0.5), STPoint(1.4, 0.5));

  /* s
   * ^
   * |
   * |
   * |---------*-----------------
   * |       *                  |
   * |     *   *                |
   * |   *   *                  |
   * | *   *                    |
   * *   *                      |
   * | *                        |
   * *--------------------------*------------------------------------> t
   *
   */
  EXPECT_EQ(point_pairs.size(), 5);

  StBoundary st_boundary;
  st_boundary.RemoveRedundantPoints(&point_pairs);

  EXPECT_EQ(point_pairs.size(), 2);
  EXPECT_DOUBLE_EQ(point_pairs[0].first.s(), 0.0);
  EXPECT_DOUBLE_EQ(point_pairs[0].first.t(), 0.0);
  EXPECT_DOUBLE_EQ(point_pairs[0].second.s(), 1.0);
  EXPECT_DOUBLE_EQ(point_pairs[0].second.t(), 0.0);

  EXPECT_DOUBLE_EQ(point_pairs[1].first.s(), 0.4);
  EXPECT_DOUBLE_EQ(point_pairs[1].first.t(), 0.5);
  EXPECT_DOUBLE_EQ(point_pairs[1].second.s(), 1.4);
  EXPECT_DOUBLE_EQ(point_pairs[1].second.t(), 0.5);
}

}  // namespace planning
}  // namespace apollo
