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
 * @file : spline_1d_constraint.cc
 * @brief: wrapp up solver constraint interface with direct methods and preset
 *methods
 **/

#include "modules/planning/math/smoothing_spline/spline_1d_constraint.h"

#include <limits>

#include "modules/common/log.h"

namespace apollo {
namespace planning {

using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;

Spline1dConstraint::Spline1dConstraint(const Spline1d& pss)
    : Spline1dConstraint(pss.x_knots(), pss.spline_order()) {}

Spline1dConstraint::Spline1dConstraint(const std::vector<double>& x_knots,
                                       const uint32_t             spline_order)
    : x_knots_(x_knots)
    , spline_order_(spline_order) {
  inequality_constraint_.SetIsEquality(false);
  equality_constraint_.SetIsEquality(true);
}

// 从上往下看AddConstraint()，首先是spline_constraint->AddInequalityConstraint()。
// 其中又内部调用了AffineConstraint::AddConstraint()。AffineConstraint::AddConstraint()是最底层的函数，简明清晰。
// 其他函数之所以复杂，正在于要构造这个函数的2个输入参数：constraint_matrix 和 constraint_boundary。
bool Spline1dConstraint::AddInequalityConstraint(const Eigen::MatrixXd& constraint_matrix,
                                                 const Eigen::MatrixXd& constraint_boundary) {
  return inequality_constraint_.AddConstraint(constraint_matrix, constraint_boundary);
}

bool Spline1dConstraint::AddEqualityConstraint(const Eigen::MatrixXd& constraint_matrix,
                                               const Eigen::MatrixXd& constraint_boundary) {
  return equality_constraint_.AddConstraint(constraint_matrix, constraint_boundary);
}

// 添加lower_bound <= f(x_coord) <= upper_bound取值范围约束
bool Spline1dConstraint::AddBoundary(const std::vector<double>& x_coord,
                                     const std::vector<double>& lower_bound,
                                     const std::vector<double>& upper_bound) {
  // 检验lower_bound和upper_bound中的元素值是否是有效值
  // 将一一对应的x_coord -> lower_bound存入 filtered_lower_bound_x 和 filtered_lower_bound
  // 将一一对应的x_coord -> upper_bound存入 filtered_upper_bound_x 和 filtered_upper_bound
  std::vector<double> filtered_lower_bound;
  std::vector<double> filtered_upper_bound;
  std::vector<double> filtered_lower_bound_x;
  std::vector<double> filtered_upper_bound_x;

  FilterConstraints(x_coord, lower_bound, upper_bound, &filtered_lower_bound_x,
                    &filtered_lower_bound, &filtered_upper_bound_x, &filtered_upper_bound);

  // emplace affine constraints
  const uint32_t num_params = spline_order_ + 1;
  // inequality_constraint矩阵中绝大多数元素都是0，
  //只有x_coord坐落的一段curve对应元素有值
  Eigen::MatrixXd inequality_constraint =
      //                                    //行              //列
      Eigen::MatrixXd::Zero(filtered_upper_bound.size() + filtered_lower_bound.size(),
                            (x_knots_.size() - 1) * num_params);

  /*
   * std::vector<double> x_coord     = {0,   5,   10,  20};
   * std::vector<double> lower_bound = {0,   0,   0,   0};
   * std::vector<double> upper_bound = {3.5, 3.5, 3.5, 3.5};
   *
   *    ^ d
   *    |
   *    | * knot point
   *    | @ constraint point
   *    |------------------------------------------ upper_bound
   *    |
   *    |0           5            10              20
   *    |*--@--------*-@--------@-*--------@------*----------> s
   *    |
   *    |
   *    |------------------------------------------ lower_bound
   *
   * Ax >= b
   * inequality_constraint:
   * | 1, s0, s0^2, s0^3, s0^4, s0^5 |         >=  |l0|
   * | 1, s1, s1^2, s1^3, s1^4, s1^5 |         >=  |l1|
   * | 1, s2, s2^2, s2^3, s2^4, s2^5 |   |a0|  >=  |l2|
   * | 1, s3, s3^2, s3^3, s3^4, s3^5 |   |a1|  >=  |l3|
   * | 1, s4, s4^2, s4^3, s4^4, s4^5 | * |a2|  >=  |l4|
   * | 1, s5, s5^2, s5^3, s5^4, s5^5 |   |a3|  >=  |..|
   * | 1, s6, s6^2, s6^3, s6^4, s6^5 |   |a4|  >=  |..|
   * | 1, s7, s7^2, s7^3, s7^4, s7^5 |   |a5|  >=  |..|
   * | 1, s8, s8^2, s8^3, s8^4, s8^5 |         >=  |..|
   * | 1, s9, s9^2, s9^3, s9^4, s9^5 |         >=  |ln|
   */

  Eigen::MatrixXd inequality_boundary =
      Eigen::MatrixXd::Zero(filtered_upper_bound.size() + filtered_lower_bound.size(), 1);

  // 由 AddPointConstraintInRange() 推来，其实就是已知val=f(x0)，range>0，
  // x_coord=filtered_lower_bound_x = filtered_upper_bound_x，
  // filtered_lower_bound=lower_bound = val-range，
  // filtered_upper_bound=upper_bound = val+range，
  // 由val-range<val<val+range ==> f(filtered_lower_bound_x) >= filtered_lower_bound
  // 且 f(filtered_upper_bound_x) <= filtered_upper_bound

  // 这个for循环对应f(filtered_lower_bound_x) >= filtered_lower_bound的情况，
  // 恰好符合QP中 Ax>=b 的形式
  for (uint32_t i = 0; i < filtered_lower_bound.size(); ++i) {
    // FindIndex()返回x_knots_中第一个大于入参的元素的前一个位置，
    //即入参位于index～index+1之间
    uint32_t index = FindIndex(filtered_lower_bound_x[i]);

    // corrected_x 是每一段的相对起点的x坐标
    const double corrected_x = filtered_lower_bound_x[i] - x_knots_[index];
    double       coef        = 1.0;

    // j=0~num_params，依次是  corrected_x 的0~num_params次项
    //计算[1, s, s^2, s^3, s^4, s^5] * [a0, a1, a2, a3, a4, a5].T，
    //[a, b, c, d, e, f] 是5次curve的系数
    for (uint32_t j = 0; j < num_params; ++j) {
      //对特定行、特定curve的参数列更新，因为给定一个s，其根据定义域只对应一段curve
      inequality_constraint(i, j + index * num_params) = coef;
      coef *= corrected_x;
      /*
       * 1
       * x
       * x^2
       * x^3
       * x^4
       * x^5
       */
    }
    inequality_boundary(i, 0) = filtered_lower_bound[i];
  }

  //这个for循环对应f(filtered_upper_bound_x) <= filtered_upper_bound的情况，
  //为符合QP中Ax>=b的形式,将Ax<=b ==> -Ax>=-b，故coef赋初值-1，与上部分比，全为负数
  for (uint32_t i = 0; i < filtered_upper_bound.size(); ++i) {
    uint32_t     index       = FindIndex(filtered_upper_bound_x[i]);
    const double corrected_x = filtered_upper_bound_x[i] - x_knots_[index];
    double       coef        = -1.0;
    //-f(x) = [-1, -x, -x^2, -x^3, -x^4, -x^5] * [a, b, c, d, e, f].T
    // >= -upper_bound ==>
    // f(x) = [1, x, x^2, x^3, x^4, x^5] * [a, b, c, d, e, f].T
    // <= upper_bound
    for (uint32_t j = 0; j < num_params; ++j) {
      // 注意行列坐标
      inequality_constraint(i + filtered_lower_bound.size(), j + index * num_params) = coef;
      coef *= corrected_x;
    }
    inequality_boundary(i + filtered_lower_bound.size(), 0) = -filtered_upper_bound[i];
  }
  // Ax >= b
  // A: inequality_constraint
  // b: inequality_boundary
  return inequality_constraint_.AddConstraint(inequality_constraint, inequality_boundary);
}

bool Spline1dConstraint::AddDerivativeBoundary(const std::vector<double>& x_coord,
                                               const std::vector<double>& lower_bound,
                                               const std::vector<double>& upper_bound) {
  std::vector<double> filtered_lower_bound;
  std::vector<double> filtered_upper_bound;
  std::vector<double> filtered_lower_bound_x;
  std::vector<double> filtered_upper_bound_x;

  if (x_knots_.size() < 2) { return false; }

  if (!FilterConstraints(x_coord, lower_bound, upper_bound, &filtered_lower_bound_x,
                         &filtered_lower_bound, &filtered_upper_bound_x, &filtered_upper_bound)) {
    return false;
  }

  // emplace affine constraints
  const uint32_t  num_params = spline_order_ + 1;
  Eigen::MatrixXd inequality_constraint =
      Eigen::MatrixXd::Zero(filtered_upper_bound.size() + filtered_lower_bound.size(),
                            (x_knots_.size() - 1) * num_params);
  Eigen::MatrixXd inequality_boundary =
      Eigen::MatrixXd::Zero(filtered_upper_bound.size() + filtered_lower_bound.size(), 1);

  for (uint32_t i = 0; i < filtered_lower_bound.size(); ++i) {
    uint32_t     index       = FindIndex(filtered_lower_bound_x[i]);
    const double corrected_x = filtered_lower_bound_x[i] - x_knots_[index];
    double       coef        = 1.0;
    // 计算[0， 1, 2s, 3s^2, 4s^3, 5s^4] * [a0, a1, a2, a3, a4, a5].T，前者是f'(s)，
    // 后者是5次curve的系数
    for (uint32_t j = 1; j < num_params; ++j) {
      inequality_constraint(i, j + index * num_params) = coef * j;
      coef *= corrected_x;
    }
    inequality_boundary(i, 0) = filtered_lower_bound[i];
  }

  for (uint32_t i = 0; i < filtered_upper_bound.size(); ++i) {
    uint32_t     index       = FindIndex(filtered_upper_bound_x[i]);
    const double corrected_x = filtered_upper_bound_x[i] - x_knots_[index];
    double       coef        = -1.0;
    for (uint32_t j = 1; j < num_params; ++j) {
      inequality_constraint(i + filtered_lower_bound.size(), j + index * num_params) = coef * j;
      coef *= corrected_x;
    }
    inequality_boundary(i + filtered_lower_bound.size(), 0) = -filtered_upper_bound[i];
  }
  return inequality_constraint_.AddConstraint(inequality_constraint, inequality_boundary);
}

bool Spline1dConstraint::AddSecondDerivativeBoundary(const std::vector<double>& x_coord,
                                                     const std::vector<double>& lower_bound,
                                                     const std::vector<double>& upper_bound) {
  std::vector<double> filtered_lower_bound;
  std::vector<double> filtered_upper_bound;
  std::vector<double> filtered_lower_bound_x;
  std::vector<double> filtered_upper_bound_x;

  if (x_knots_.size() < 2) { return false; }

  if (!FilterConstraints(x_coord, lower_bound, upper_bound, &filtered_lower_bound_x,
                         &filtered_lower_bound, &filtered_upper_bound_x, &filtered_upper_bound)) {
    return false;
  }

  // emplace affine constraints
  const uint32_t  num_params = spline_order_ + 1;
  Eigen::MatrixXd inequality_constraint =
      Eigen::MatrixXd::Zero(filtered_upper_bound.size() + filtered_lower_bound.size(),
                            (x_knots_.size() - 1) * num_params);
  Eigen::MatrixXd inequality_boundary =
      Eigen::MatrixXd::Zero(filtered_upper_bound.size() + filtered_lower_bound.size(), 1);

  for (uint32_t i = 0; i < filtered_lower_bound.size(); ++i) {
    uint32_t     index       = FindIndex(filtered_lower_bound_x[i]);
    const double corrected_x = filtered_lower_bound_x[i] - x_knots_[index];
    double       coef        = 1.0;
    //计算[0, 0, 2, 6s, 12s^2, 20s^3] * [a0, a1, a2, a3, a4, a5].T，前者是f''(s),
    //后者是5次curve的系数
    for (uint32_t j = 2; j < num_params; ++j) {
      inequality_constraint(i, j + index * num_params) = coef * j * (j - 1);
      coef *= corrected_x;
    }
    inequality_boundary(i, 0) = filtered_lower_bound[i];
  }

  for (uint32_t i = 0; i < filtered_upper_bound.size(); ++i) {
    uint32_t     index       = FindIndex(filtered_upper_bound_x[i]);
    const double corrected_x = filtered_upper_bound_x[i] - x_knots_[index];
    double       coef        = -1.0;
    for (uint32_t j = 2; j < num_params; ++j) {
      inequality_constraint(i + filtered_lower_bound.size(), j + index * num_params) =
          coef * j * (j - 1);
      coef *= corrected_x;
    }
    inequality_boundary(i + filtered_lower_bound.size(), 0) = -filtered_upper_bound[i];
  }
  return inequality_constraint_.AddConstraint(inequality_constraint, inequality_boundary);
}

bool Spline1dConstraint::AddThirdDerivativeBoundary(const std::vector<double>& x_coord,
                                                    const std::vector<double>& lower_bound,
                                                    const std::vector<double>& upper_bound) {
  std::vector<double> filtered_lower_bound;
  std::vector<double> filtered_upper_bound;
  std::vector<double> filtered_lower_bound_x;
  std::vector<double> filtered_upper_bound_x;

  if (!FilterConstraints(x_coord, lower_bound, upper_bound, &filtered_lower_bound_x,
                         &filtered_lower_bound, &filtered_upper_bound_x, &filtered_upper_bound)) {
    AERROR << "Fail to filter constraints.";
    return false;
  }

  if (x_knots_.size() < 2) {
    AERROR << "x_konts size cannot be < 2.";
    return false;
  }

  // emplace affine constraints
  const uint32_t  num_params = spline_order_ + 1;
  Eigen::MatrixXd inequality_constraint =
      Eigen::MatrixXd::Zero(filtered_upper_bound.size() + filtered_lower_bound.size(),
                            (x_knots_.size() - 1) * num_params);
  Eigen::MatrixXd inequality_boundary =
      Eigen::MatrixXd::Zero(filtered_upper_bound.size() + filtered_lower_bound.size(), 1);

  for (uint32_t i = 0; i < filtered_lower_bound.size(); ++i) {
    uint32_t     index       = FindIndex(filtered_lower_bound_x[i]);
    const double corrected_x = filtered_lower_bound_x[i] - x_knots_[index];
    double       coef        = 1.0;
    for (uint32_t j = 3; j < num_params; ++j) {
      inequality_constraint(i, j + index * num_params) = coef * j * (j - 1) * (j - 2);
      coef *= corrected_x;
    }
    inequality_boundary(i, 0) = filtered_lower_bound[i];
  }

  for (uint32_t i = 0; i < filtered_upper_bound.size(); ++i) {
    uint32_t     index       = FindIndex(filtered_upper_bound_x[i]);
    const double corrected_x = filtered_upper_bound_x[i] - x_knots_[index];
    double       coef        = -1.0;
    for (uint32_t j = 3; j < num_params; ++j) {
      inequality_constraint(i + filtered_lower_bound.size(), j + index * num_params) =
          coef * j * (j - 1) * (j - 2);
      coef *= corrected_x;
    }
    inequality_boundary(i + filtered_lower_bound.size(), 0) = -filtered_upper_bound[i];
  }
  return inequality_constraint_.AddConstraint(inequality_constraint, inequality_boundary);
}

/*
 *
 * |1   x   x^2   x^3  x^4  x^5| * |c0  c1  c2  c3  c4  c5|^T = fx
 *
 *
 */
bool Spline1dConstraint::AddPointConstraint(const double x, const double fx) {
  uint32_t            index = FindIndex(x);
  std::vector<double> power_x;
  const uint32_t      num_params = spline_order_ + 1;
  GeneratePowerX(x - x_knots_[index], num_params, &power_x);
  Eigen::MatrixXd equality_constraint =
      Eigen::MatrixXd::Zero(1, (x_knots_.size() - 1) * num_params);
  uint32_t index_offset = index * num_params;
  for (uint32_t i = 0; i < num_params; ++i) {
    equality_constraint(0, index_offset + i) = power_x[i];
  }
  Eigen::MatrixXd equality_boundary(1, 1);
  equality_boundary(0, 0) = fx;
  return AddEqualityConstraint(equality_constraint, equality_boundary);
}

bool Spline1dConstraint::AddConstraintInRange(AddConstraintInRangeFunc func,
                                              const double             x,
                                              const double             val,
                                              const double             range) {
  if (range < 0.0) { return false; }
  std::vector<double> x_vec;
  x_vec.push_back(x);

  std::vector<double> lower_bound;
  std::vector<double> upper_bound;
  lower_bound.push_back(val - range);
  upper_bound.push_back(val + range);
  return func(x_vec, lower_bound, upper_bound);
}
/*
 * 添加0阶不等式约束，即考察目标函数f(x)与 [fx-range, fx+range] 的不等式关系。
 * 间接调用的Spline1dConstraint::AddBoundary()是理解constraint的核心，
 * 它就是在构造QP形式中Ax>=b约束中的A和b
 *
 *
 * |1   x   x^2   x^3     x^4       x^5| * |c0  c1  c2  c3  c4  c5|^T = fx
 * |0   1   2*x   3*x^2   4*x^3   5*x^4| * |c0  c1  c2  c3  c4  c5|^T = fx
 *
 *
 */
bool Spline1dConstraint::AddPointConstraintInRange(const double x,
                                                   const double fx,
                                                   const double range) {
  return AddConstraintInRange(std::bind(&Spline1dConstraint::AddBoundary, this, _1, _2, _3), x, fx,
                              range);
}

bool Spline1dConstraint::AddPointDerivativeConstraint(const double x, const double dfx) {
  uint32_t            index = FindIndex(x);
  std::vector<double> power_x;
  const uint32_t      num_params = spline_order_ + 1;
  GeneratePowerX(x - x_knots_[index], num_params, &power_x);
  Eigen::MatrixXd equality_constraint =
      Eigen::MatrixXd::Zero(1, (x_knots_.size() - 1) * num_params);
  uint32_t index_offset = index * num_params;
  for (uint32_t i = 1; i < num_params; ++i) {
    //该行矩阵中对应该spline seg是[0, 1, 2*x, 3*x^2, 4*x^3, 5*x^4]
    //注意equality_constraint的长度并不是6，此处只分析某一段的长度是6
    //注意equality_constraint被初始化为元素全为0的行矩阵
    equality_constraint(0, index_offset + i) = power_x[i - 1] * i;
  }
  Eigen::MatrixXd equality_boundary(1, 1);
  equality_boundary(0, 0) = dfx;
  //添加在某一点的一阶导数相等约束，对5次多项式曲线，如下
  //[0, 1, 2*x, 3*x^2, 4*x^3, 5*x^4] * [a, b, c, d, e, f].T = dfx
  //中间是5次curve的系数，后者是该点处的一阶导数值dfx
  return AddEqualityConstraint(equality_constraint, equality_boundary);
}

// 添加1阶不等式约束（angle），即考察目标函数f'(x)与 [dfx-range, dfx+range] 的不等式关系。
// 间接调用的Spline1dConstraint::AddDerivativeBoundary()与Spline1dConstraint::AddBoundary()如出一辙。
// 值得注意的是，此处的A是f(x)的一阶导。
bool Spline1dConstraint::AddPointDerivativeConstraintInRange(const double x,
                                                             const double dfx,
                                                             const double range) {
  return AddConstraintInRange(
      std::bind(&Spline1dConstraint::AddDerivativeBoundary, this, _1, _2, _3), x, dfx, range);
}

bool Spline1dConstraint::AddPointSecondDerivativeConstraint(const double x, const double ddfx) {
  uint32_t            index = FindIndex(x);
  std::vector<double> power_x;
  const uint32_t      num_params = spline_order_ + 1;
  GeneratePowerX(x - x_knots_[index], num_params, &power_x);
  Eigen::MatrixXd equality_constraint =
      Eigen::MatrixXd::Zero(1, (x_knots_.size() - 1) * num_params);
  uint32_t index_offset = index * num_params;
  for (uint32_t i = 2; i < num_params; ++i) {
    equality_constraint(0, index_offset + i) = power_x[i - 2] * i * (i - 1);
  }
  Eigen::MatrixXd equality_boundary(1, 1);
  equality_boundary(0, 0) = ddfx;
  return AddEqualityConstraint(equality_constraint, equality_boundary);
}

// 添加2阶不等式约束（kappa）。间接调用的Spline1dConstraint::AddSecondDerivativeBoundary()的A是f(x)的二阶导。
bool Spline1dConstraint::AddPointSecondDerivativeConstraintInRange(const double x,
                                                                   const double ddfx,
                                                                   const double range) {
  return AddConstraintInRange(
      std::bind(&Spline1dConstraint::AddSecondDerivativeBoundary, this, _1, _2, _3), x, ddfx,
      range);
}

bool Spline1dConstraint::AddPointThirdDerivativeConstraint(const double x, const double dddfx) {
  uint32_t            index = FindIndex(x);
  std::vector<double> power_x;
  const uint32_t      num_params = spline_order_ + 1;
  GeneratePowerX(x - x_knots_[index], num_params, &power_x);
  Eigen::MatrixXd equality_constraint =
      Eigen::MatrixXd::Zero(1, (x_knots_.size() - 1) * num_params);
  uint32_t index_offset = index * num_params;
  for (uint32_t i = 3; i < num_params; ++i) {
    equality_constraint(0, index_offset + i) = power_x[i - 3] * i * (i - 1) * (i - 2);
  }
  Eigen::MatrixXd equality_boundary(1, 1);
  equality_boundary(0, 0) = dddfx;
  return AddEqualityConstraint(equality_constraint, equality_boundary);
}

bool Spline1dConstraint::AddPointThirdDerivativeConstraintInRange(const double x,
                                                                  const double dddfx,
                                                                  const double range) {
  return AddConstraintInRange(
      std::bind(&Spline1dConstraint::AddSecondDerivativeBoundary, this, _1, _2, _3), x, dddfx,
      range);
}

bool Spline1dConstraint::AddSmoothConstraint() {
  if (x_knots_.size() < 3) { return false; }
  const uint32_t  num_params = spline_order_ + 1;
  Eigen::MatrixXd equality_constraint =
      Eigen::MatrixXd::Zero(x_knots_.size() - 2, (x_knots_.size() - 1) * num_params);
  Eigen::MatrixXd equality_boundary = Eigen::MatrixXd::Zero(x_knots_.size() - 2, 1);

  for (uint32_t i = 0; i < x_knots_.size() - 2; ++i) {
    double       left_coef  = 1.0;
    double       right_coef = -1.0;
    const double x_left     = x_knots_[i + 1] - x_knots_[i];
    const double x_right    = 0.0;
    for (uint32_t j = 0; j < num_params; ++j) {
      equality_constraint(i, num_params * i + j)       = left_coef;
      equality_constraint(i, num_params * (i + 1) + j) = right_coef;
      left_coef *= x_left;
      right_coef *= x_right;
    }
  }
  return equality_constraint_.AddConstraint(equality_constraint, equality_boundary);
}

bool Spline1dConstraint::AddDerivativeSmoothConstraint() {
  if (x_knots_.size() < 3) { return false; }

  const uint32_t  n_constraint = (x_knots_.size() - 2) * 2;
  const uint32_t  num_params   = spline_order_ + 1;
  Eigen::MatrixXd equality_constraint =
      Eigen::MatrixXd::Zero(n_constraint, (x_knots_.size() - 1) * num_params);
  Eigen::MatrixXd equality_boundary = Eigen::MatrixXd::Zero(n_constraint, 1);

  for (uint32_t i = 0; i < n_constraint; i += 2) {
    double       left_coef   = 1.0;
    double       right_coef  = -1.0;
    double       left_dcoef  = 1.0;
    double       right_dcoef = -1.0;
    const double x_left      = x_knots_[i / 2 + 1] - x_knots_[i / 2];
    const double x_right     = 0.0;
    for (uint32_t j = 0; j < num_params; ++j) {
      equality_constraint(i, num_params * (i / 2) + j)       = left_coef;
      equality_constraint(i, num_params * ((i / 2) + 1) + j) = right_coef;
      if (j >= 1) {
        equality_constraint(i + 1, num_params * (i / 2) + j)       = left_dcoef * j;
        equality_constraint(i + 1, num_params * ((i / 2) + 1) + j) = right_dcoef * j;
        left_dcoef                                                 = left_coef;
        right_dcoef                                                = right_coef;
      }
      left_coef *= x_left;
      right_coef *= x_right;
    }
  }
  return equality_constraint_.AddConstraint(equality_constraint, equality_boundary);
}

bool Spline1dConstraint::AddSecondDerivativeSmoothConstraint() {
  if (x_knots_.size() < 3) { return false; }

  const uint32_t  n_constraint = (x_knots_.size() - 2) * 3;
  const uint32_t  num_params   = spline_order_ + 1;
  Eigen::MatrixXd equality_constraint =
      Eigen::MatrixXd::Zero(n_constraint, (x_knots_.size() - 1) * num_params);
  Eigen::MatrixXd equality_boundary = Eigen::MatrixXd::Zero(n_constraint, 1);

  for (uint32_t i = 0; i < n_constraint; i += 3) {
    double left_coef    = 1.0;
    double right_coef   = -1.0;
    double left_dcoef   = 1.0;
    double right_dcoef  = -1.0;
    double left_ddcoef  = 1.0;
    double right_ddcoef = -1.0;

    const double x_left  = x_knots_[i / 3 + 1] - x_knots_[i / 3];
    const double x_right = 0.0;
    for (uint32_t j = 0; j < num_params; ++j) {
      equality_constraint(i, num_params * (i / 3) + j)     = left_coef;
      equality_constraint(i, num_params * (i / 3 + 1) + j) = right_coef;

      if (j >= 2) {
        equality_constraint(i + 2, num_params * i / 3 + j)       = left_ddcoef * j * (j - 1);
        equality_constraint(i + 2, num_params * (i / 3 + 1) + j) = right_ddcoef * j * (j - 1);
        left_ddcoef                                              = left_dcoef;
        right_ddcoef                                             = right_dcoef;
      }

      if (j >= 1) {
        equality_constraint(i + 1, num_params * (i / 3) + j)     = left_dcoef * j;
        equality_constraint(i + 1, num_params * (i / 3 + 1) + j) = right_dcoef * j;
        left_dcoef                                               = left_coef;
        right_dcoef                                              = right_coef;
      }
      left_coef *= x_left;
      right_coef *= x_right;
    }
  }
  return equality_constraint_.AddConstraint(equality_constraint, equality_boundary);
}

// 添加joint points处0~3阶导连续 等式约束，即
// fl(xl)=fr(xr); fl'(xl)=fr'(xr); fl''(xl)=fr''(xr); fl'''(xl)=fr'''(xr)，
// 转换为QP的形式，即
// fl(xl)-fr(xr)=0; fl'(xl)-fr'(xr)=0; fl''(xl)-fr''(xr)=0; fl'''(xl)-fr'''(xr)=0
// fl()和fr()是相邻的2段curve polynomial，
// xl和xr是同一个点在不同的段不同的起点下的相对坐标
bool Spline1dConstraint::AddThirdDerivativeSmoothConstraint() {
  if (x_knots_.size() < 3) { return false; }

  const uint32_t  n_constraint = (x_knots_.size() - 2) * 4;
  const uint32_t  num_params   = spline_order_ + 1;
  Eigen::MatrixXd equality_constraint =
      Eigen::MatrixXd::Zero(n_constraint, (x_knots_.size() - 1) * num_params);
  Eigen::MatrixXd equality_boundary = Eigen::MatrixXd::Zero(n_constraint, 1);
  // 这里的2层for循环很不直观，其实是把不同阶导数、不同的joint point、不同的系数
  // 杂糅在一起计算了。虽高效，但太复杂抽象了，最终计算后应该是下面的形式
  // 每行左边6项是joint point左侧的curve参数，右边6项是joint point右侧的curve参数
  // 0~3行是1个joint point的0~3阶导，每4行表示一个点
  // | 1 xl xl^2 xl^3  xl^4   xl^5   -1 -xr -xr^2 -xr^3  -xr^4   -xr^5   |   | al | = | 0 |
  // | 0 1  2xl  3xl^2 4xl^3  5xl^4   0 -1  -2xr  -3xr^2 -4xr^3  -5xr^4  |   | bl | = | 0 |
  // | 0 0  2    6xl   12xl^2 20xl^3  0  0  -2    -6xr   -12xr^2 -20xr^3 | * | cl | = | 0 |
  // | 0 0  0    6     24xl   60xl^2  0  0   0    -6     -24xr   -60xr^2 |   | dl | = | 0 |
  // | . .  .    .     .      .       .  .   .     .      .       .      |   | el | = | 0 |
  // | . .  .    .     .      .       .  .   .     .      .       .      |   | fl | = | 0 |
  //                                                                         | ar | = | 0 |
  //                                                                         | br | = | 0 |
  //                                                                         | cr | = | 0 |
  //                                                                         | dr | = | 0 |
  //                                                                         | er | = | 0 |
  //                                                                         | fr | = | 0 |

  // i循环x_knots_.size())-2次
  for (uint32_t i = 0; i < n_constraint; i += 4) {
    double left_coef     = 1.0;
    double right_coef    = -1.0;
    double left_dcoef    = 1.0;
    double right_dcoef   = -1.0;
    double left_ddcoef   = 1.0;
    double right_ddcoef  = -1.0;
    double left_dddcoef  = 1.0;
    double right_dddcoef = -1.0;

    //第n段的最后一个点相对于第n段的第一点的坐标
    const double x_left = x_knots_[i / 4 + 1] - x_knots_[i / 4];
    //第n段的最后一个点相对于第n+1段的第一点的坐标
    const double x_right = 0.0;
    for (uint32_t j = 0; j < num_params; ++j) {
      // j循环6次
      equality_constraint(i, num_params * i / 4 + j)       = left_coef;
      equality_constraint(i, num_params * (i / 4 + 1) + j) = right_coef;

      if (j >= 3) {
        equality_constraint(i + 3, num_params * i / 4 + j) = left_dddcoef * j * (j - 1) * (j - 2);
        equality_constraint(i + 3, num_params * (i / 4 + 1) + j) =
            right_dddcoef * j * (j - 1) * (j - 2);
        left_dddcoef  = left_ddcoef;
        right_dddcoef = right_ddcoef;
      }

      if (j >= 2) {
        equality_constraint(i + 2, num_params * i / 4 + j)       = left_ddcoef * j * (j - 1);
        equality_constraint(i + 2, num_params * (i / 4 + 1) + j) = right_ddcoef * j * (j - 1);
        left_ddcoef                                              = left_dcoef;
        right_ddcoef                                             = right_dcoef;
      }

      if (j >= 1) {
        equality_constraint(i + 1, num_params * i / 4 + j)       = left_dcoef * j;
        equality_constraint(i + 1, num_params * (i / 4 + 1) + j) = right_dcoef * j;
        left_dcoef                                               = left_coef;
        right_dcoef                                              = right_coef;
      }

      left_coef *= x_left;
      right_coef *= x_right;
    }
  }
  return equality_constraint_.AddConstraint(equality_constraint, equality_boundary);
}

// 单调性，表示下个点的 s值 要比它的前一个s值大，表示，车是处于前进状态
bool Spline1dConstraint::AddMonotoneInequalityConstraint(const std::vector<double>& x_coord) {
  if (x_coord.size() < 2) {
    // Skip because NO inequality constraint is needed.
    return false;
  }

  const uint32_t num_params = spline_order_ + 1;
  //注意维度，有N个x，就有N-1个不等式f(xn-1) <= f(xn)
  Eigen::MatrixXd inequality_constraint =
      Eigen::MatrixXd::Zero(x_coord.size() - 1, (x_knots_.size() - 1) * num_params);

  // size()-1 行 1 列 个 0， 表示连续两个点的 s 差值大于0 (cur_coef[j] - prev_coef[j] >= 0)
  Eigen::MatrixXd inequality_boundary = Eigen::MatrixXd::Zero(x_coord.size() - 1, 1);

  uint32_t            prev_spline_index = FindIndex(x_coord[0]);
  double              prev_rel_x        = x_coord[0] - x_knots_[prev_spline_index];
  std::vector<double> prev_coef;
  // s = s0 + s1 + s2^2 + s3^3 + s4^4 + s5^5
  GeneratePowerX(prev_rel_x, num_params, &prev_coef);
  for (uint32_t i = 1; i < x_coord.size(); ++i) {
    uint32_t            cur_spline_index = FindIndex(x_coord[i]);
    double              cur_rel_x        = x_coord[i] - x_knots_[cur_spline_index];
    std::vector<double> cur_coef;

    GeneratePowerX(cur_rel_x, num_params, &cur_coef);
    // if constraint on the same spline
    if (cur_spline_index == prev_spline_index) {
      for (uint32_t j = 0; j < cur_coef.size(); ++j) {
        //在同一段spline，则多项式系数是相同的。
        //记cur_rel_x-prev_rel_x=d, inequality_constraint对应该spline seg是
        //[1-1, d, d^2, d^3, d^4, d^5] * [a, b, c, d, e, f].T >= 0
        inequality_constraint(i - 1, cur_spline_index * num_params + j) =
            cur_coef[j] - prev_coef[j];
      }
    } else {
      //[-1, -px, -px^2, -px^3, -px^4, -px^5, 1, cx, cx^2, cx^3, cx^4, cx^5]
      // * [pa, pb, pc, pd, pe, pf, ca, cb, cc, cd, ce, cf].T >= 0
      // if not on the same spline
      for (uint32_t j = 0; j < cur_coef.size(); ++j) {
        inequality_constraint(i - 1, prev_spline_index * num_params + j) = -prev_coef[j];
        inequality_constraint(i - 1, cur_spline_index * num_params + j)  = cur_coef[j];
      }
    }
    prev_coef         = cur_coef;
    prev_spline_index = cur_spline_index;
  }

  return inequality_constraint_.AddConstraint(inequality_constraint, inequality_boundary);
}

bool Spline1dConstraint::AddMonotoneInequalityConstraintAtKnots() {
  return AddMonotoneInequalityConstraint(x_knots_);
}

const AffineConstraint& Spline1dConstraint::inequality_constraint() const {
  return inequality_constraint_;
}

const AffineConstraint& Spline1dConstraint::equality_constraint() const {
  return equality_constraint_;
}

uint32_t Spline1dConstraint::FindIndex(const double x) const {
  auto upper_bound = std::upper_bound(x_knots_.begin() + 1, x_knots_.end(), x);
  return std::min(static_cast<uint32_t>(x_knots_.size() - 1),
                  static_cast<uint32_t>(upper_bound - x_knots_.begin())) -
         1;
}

bool Spline1dConstraint::FilterConstraints(const std::vector<double>& x_coord,
                                           const std::vector<double>& lower_bound,
                                           const std::vector<double>& upper_bound,
                                           std::vector<double>* const filtered_lower_bound_x,
                                           std::vector<double>* const filtered_lower_bound,
                                           std::vector<double>* const filtered_upper_bound_x,
                                           std::vector<double>* const filtered_upper_bound) {
  filtered_lower_bound->clear();
  filtered_upper_bound->clear();
  filtered_lower_bound_x->clear();
  filtered_upper_bound_x->clear();

  const double inf = std::numeric_limits<double>::infinity();

  filtered_lower_bound->reserve(lower_bound.size());
  filtered_lower_bound_x->reserve(lower_bound.size());

  filtered_upper_bound->reserve(upper_bound.size());
  filtered_upper_bound_x->reserve(upper_bound.size());

  for (uint32_t i = 0; i < lower_bound.size(); ++i) {
    if (std::isnan(lower_bound[i]) || lower_bound[i] == inf) { return false; }
    if (lower_bound[i] < inf && lower_bound[i] > -inf) {
      filtered_lower_bound->emplace_back(lower_bound[i]);
      filtered_lower_bound_x->emplace_back(x_coord[i]);
    }
  }

  for (uint32_t i = 0; i < upper_bound.size(); ++i) {
    if (std::isnan(upper_bound[i]) || upper_bound[i] == -inf) { return false; }
    if (upper_bound[i] < inf && upper_bound[i] > -inf) {
      filtered_upper_bound->emplace_back(upper_bound[i]);
      filtered_upper_bound_x->emplace_back(x_coord[i]);
    }
  }
  return true;
}

void Spline1dConstraint::GeneratePowerX(const double               x,
                                        const uint32_t             order,
                                        std::vector<double>* const power_x) const {
  double cur_x = 1.0;
  for (uint32_t i = 0; i < order; ++i) {
    power_x->push_back(cur_x);
    cur_x *= x;
  }
}

}  // namespace planning
}  // namespace apollo
