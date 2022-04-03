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
 * @brief Defines the Vec2d class.
 */

#ifndef MODULES_COMMON_MATH_VEC2D_H_
#define MODULES_COMMON_MATH_VEC2D_H_

#include <cmath>
#include <string>

/**
 * @namespace apollo::common::math
 * @brief apollo::common::math
 */
namespace apollo {
namespace common {
namespace math {

constexpr double kMathEpsilon = 1e-10;

/**
 * @class Vec2d
 *
 * @brief Implements a class of 2-dimensional vectors.
 */
// 二维向量
class Vec2d {
 public:
  //! Constructor which takes x- and y-coordinates.
  constexpr Vec2d(const double x, const double y) noexcept
      : x_(x)
      , y_(y) {}

  //! Constructor returning the zero vector.
  constexpr Vec2d() noexcept
      : Vec2d(0, 0) {}

  //! Creates a unit-vector with a given angle to the positive x semi-axis
  static Vec2d CreateUnitVec2d(const double angle);
  double       x() const { return x_; }
  double       y() const { return y_; }
  void         set_x(const double x) { x_ = x; }
  void         set_y(const double y) { y_ = y; }
  double       Length() const;
  double       LengthSquare() const;
  double       Angle() const;  // Gets the angle between the vector and the positive x semi-axis
  void         Normalize();    // Returns the unit vector that is co-linear with this vector
  double       DistanceTo(const Vec2d& other) const;        // Returns the distance to given vector
  double       DistanceSquareTo(const Vec2d& other) const;  // the squared distance to given vector
  double       CrossProd(const Vec2d& other) const;         // cross product(non-standard)
  double       InnerProd(const Vec2d& other) const;
  Vec2d        rotate(const double angle) const;
  Vec2d        operator+(const Vec2d& other) const;  // 返回一个新向量
  Vec2d        operator-(const Vec2d& other) const;  // 返回一个新向量
  Vec2d        operator*(const double ratio) const;  // 返回一个新向量
  Vec2d        operator/(const double ratio) const;  // 返回一个新向量
  Vec2d&       operator+=(const Vec2d& other);       // 返回引用(this)
  Vec2d&       operator-=(const Vec2d& other);       // 返回引用(this)
  Vec2d&       operator*=(const double ratio);       // 返回引用(this)
  Vec2d&       operator/=(const double ratio);       // 返回引用(this)
  bool         operator==(const Vec2d& other) const;
  std::string  DebugString() const;

 protected:
  double x_ = 0.0;
  double y_ = 0.0;
  // 二维向量
};

//! Multiplies the given Vec2d by a given scalar
Vec2d operator*(const double ratio, const Vec2d& vec);

}  // namespace math
}  // namespace common
}  // namespace apollo

#endif /* MODULES_COMMON_MATH_VEC2D_H_ */
