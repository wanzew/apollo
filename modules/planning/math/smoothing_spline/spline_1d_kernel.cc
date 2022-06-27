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
 * @file : piecewise_smoothing_spline_constraint.h
 * @brief: wrap up solver constraint interface with direct methods and preset
 *methods
 **/

#include "modules/planning/math/smoothing_spline/spline_1d_kernel.h"

#include <algorithm>

#include "modules/common/log.h"
#include "modules/planning/math/smoothing_spline/spline_seg_kernel.h"

namespace apollo {
namespace planning {

Spline1dKernel::Spline1dKernel(const Spline1d& spline1d)
    : Spline1dKernel(spline1d.x_knots(), spline1d.spline_order()) {}

Spline1dKernel::Spline1dKernel(const std::vector<double>& x_knots, const uint32_t spline_order)
    : x_knots_(x_knots)
    , spline_order_(spline_order) {
  // 所有curve方程的参数个数
  total_params_  = (x_knots.size() > 1 ? (x_knots.size() - 1) * (1 + spline_order_) : 0);
  kernel_matrix_ = Eigen::MatrixXd::Zero(total_params_, total_params_);
  offset_        = Eigen::MatrixXd::Zero(total_params_, 1);
}

void Spline1dKernel::AddRegularization(const double regularized_param) {
  Eigen::MatrixXd id_matrix =
      Eigen::MatrixXd::Identity(kernel_matrix_.rows(), kernel_matrix_.cols());
  kernel_matrix_ += 2.0 * id_matrix * regularized_param;
}

bool Spline1dKernel::AddKernel(const Eigen::MatrixXd& kernel,
                               const Eigen::MatrixXd& offset,
                               const double           weight) {
  if (kernel.rows() != kernel.cols() || kernel.rows() != kernel_matrix_.rows() ||
      offset.cols() != 1 || offset.rows() != offset_.rows()) {
    return false;
  }
  kernel_matrix_ += kernel * weight;
  offset_ += offset * weight;
  return true;
}

bool Spline1dKernel::AddKernel(const Eigen::MatrixXd& kernel,  //
                               const double           weight) {
  Eigen::MatrixXd offset = Eigen::MatrixXd::Zero(kernel.rows(), 1);
  return AddKernel(kernel, offset, weight);
}

Eigen::MatrixXd*       Spline1dKernel::mutable_kernel_matrix() { return &kernel_matrix_; }
Eigen::MatrixXd*       Spline1dKernel::mutable_offset() { return &offset_; }
const Eigen::MatrixXd& Spline1dKernel::kernel_matrix() const { return kernel_matrix_; }
const Eigen::MatrixXd& Spline1dKernel::offset() const { return offset_; }

// build-in kernel methods
/** 由AddNthDerivativekernelMatrix()中kernel_matrix_更新元素的方式，我们可以看出kernel_matrix_一定是下面这种形式，n
 * = spline order + 1，每一个n*n矩阵都对应着一段spline。
 * | n*n  0   0   0 |
 * |  0  n*n  0   0 |
 * |  0  0   n*n  0 |
 * |  0  0    0  n*n|
 **/
void Spline1dKernel::AddNthDerivativekernelMatrix(const uint32_t n, const double weight) {
  const uint32_t num_params = spline_order_ + 1;
  for (uint32_t i = 0; i + 1 < x_knots_.size(); ++i) {
    // weight是公式中的权重，2不理解，对于所有的cost项都×2，相当于没有影响
    // clang-format off
    Eigen::MatrixXd cur_kernel = 2 * SplineSegKernel::instance()->NthDerivativeKernel(
                                     n, 
                                     num_params, 
                                     x_knots_[i + 1] - x_knots_[i]) 
                                  * weight;
    // clang-format on
    kernel_matrix_.block(i * num_params, i * num_params, num_params, num_params) += cur_kernel;
  }
}

void Spline1dKernel::AddDerivativeKernelMatrix(const double weight) {
  AddNthDerivativekernelMatrix(1, weight);
}

void Spline1dKernel::AddSecondOrderDerivativeMatrix(const double weight) {
  AddNthDerivativekernelMatrix(2, weight);
}

void Spline1dKernel::AddThirdOrderDerivativeMatrix(const double weight) {
  AddNthDerivativekernelMatrix(3, weight);
}

void Spline1dKernel::AddNthDerivativekernelMatrixForSplineK(const uint32_t n,
                                                            const uint32_t k,
                                                            const double   weight) {
  if (k + 1 >= x_knots_.size()) {
    AERROR << "Cannot add NthDerivativeKernel for spline K because k is out of "
              "range. k = "
           << k;
    return;
  }
  const uint32_t  num_params = spline_order_ + 1;
  Eigen::MatrixXd cur_kernel = 2 *
                               SplineSegKernel::instance()->NthDerivativeKernel(
                                   n, num_params, x_knots_[k + 1] - x_knots_[k]) *
                               weight;
  kernel_matrix_.block(k * num_params, k * num_params, num_params, num_params) += cur_kernel;
}

void Spline1dKernel::AddDerivativeKernelMatrixForSplineK(const uint32_t k, const double weight) {
  AddNthDerivativekernelMatrixForSplineK(1, k, weight);
}

void Spline1dKernel::AddSecondOrderDerivativeMatrixForSplineK(const uint32_t k,
                                                              const double   weight) {
  AddNthDerivativekernelMatrixForSplineK(2, k, weight);
}

void Spline1dKernel::AddThirdOrderDerivativeMatrixForSplineK(const uint32_t k,
                                                             const double   weight) {
  AddNthDerivativekernelMatrixForSplineK(3, k, weight);
}

// 基于参考线或历史轨迹的代价
bool Spline1dKernel::AddReferenceLineKernelMatrix(const std::vector<double>& x_coord,
                                                  const std::vector<double>& ref_x,
                                                  const double               weight) {
  if (ref_x.size() != x_coord.size()) { return false; }

  /**
   * std::vector<double> x_coord     = {0,   5,   10,  20};
   * std::vector<double> lower_bound = {0,   0,   0,   0};
   * std::vector<double> upper_bound = {3.5, 3.5, 3.5, 3.5};
   *
   *    ^ d
   *    |
   *    | # knot point
   *    | * evaluated_s point
   *    |------------------------------------------ upper_bound
   *    |
   *    |0           5             10                20
   *    |*--*--*--*--#--*--*--*--*--#--*--*--*--*--*--#--> s
   *    |
   *    |
   *    |------------------------------------------ lower_bound
   *
   * Ax >= b
   * cost:
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
   **/

  const uint32_t num_params = spline_order_ + 1;
  for (uint32_t i = 0; i < x_coord.size(); ++i) {
    double cur_index = FindIndex(x_coord[i]);
    double cur_rel_x = x_coord[i] - x_knots_[cur_index];
    // update offset
    double offset_coef = -2.0 * ref_x[i] * weight;
    for (uint32_t j = 0; j < num_params; ++j) {
      offset_(j + cur_index * num_params, 0) += offset_coef;
      offset_coef *= cur_rel_x;
    }
    /**
     * 200 | 0
     * 400 | 1
     * 400 | 2
     * .......
     *     | 5
     * 200 | 0
     * 400 | 1
     * 400 | 2
     * .......
     *     | 5
     * 200 | 0
     * 400 | 1
     * 400 | 2
     * .......
     *     | 5
     * ......6*knot
     **/

    // update kernel matrix
    Eigen::MatrixXd ref_kernel(num_params, num_params);

    double              cur_x = 1.0;
    std::vector<double> power_x;
    // 循环11次 index: [0-10]
    for (uint32_t n = 0; n + 1 < 2 * num_params; ++n) {
      power_x.emplace_back(cur_x);
      cur_x *= cur_rel_x;
    }
    /* power_x:
     * [cur_x0, cur_x1, cur_x2, cur_x3, cur_x4, cur_x5, cur_x6, cur_x7, cur_x8, cur_x9, cur_x10]
     * ref_kernel:
     * | cur_x0 cur_x1 cur_x2 cur_x3 cur_x4 cur_x5 |
     * | cur_x1 cur_x2 cur_x3 cur_x4 cur_x5 cur_x6 |
     * | cur_x2 cur_x3 cur_x4 cur_x5 cur_x6 cur_x7 |
     * | cur_x3 cur_x4 cur_x5 cur_x6 cur_x7 cur_x8 |
     * | cur_x4 cur_x5 cur_x6 cur_x7 cur_x8 cur_x9 |
     * | cur_x5 cur_x6 cur_x7 cur_x8 cur_x9 cur_x10|
     *
     * e.g.
     * | 0       1       2       3       4       5 |
     * | 1       2       3       4       5       6 |
     * | 2       3       4       5       6       7 |
     * | 3       4       5       6       7       8 |
     * | 4       5       6       7       8       9 |
     * | 5       6       7       8       9       10|
     *
     *
     * 我在初次看这块代码的时候，就被csc_matrix这奇葩的矩阵构造方式折磨了好久，
     * 并且我看网上许多讲解Apollo 二次规划的文章也都没有具体到矩阵的实际形式的。
     * 我这里把代价函数P矩阵实际解构了出来，我觉得对于第一次接触这块代码的朋友，
     * 能够直观的看到这个矩阵还是会对理解Apollo的算法思想有很大帮助的
     */
    for (uint32_t r = 0; r < num_params; ++r) {
      for (uint32_t c = 0; c < num_params; ++c) {
        // 可以注意到每个元素前都乘以了2，这是为了和二次优化问题的一般形式中的1/2进行抵消的
        ref_kernel(r, c) = 2.0 * power_x[r + c];
      }
    }

    // P.block(i, j, rows, cols)          // P(i+1 : i+rows, j+1 : j+cols)
    // P.block<rows, cols>(i, j)          // P(i+1 : i+rows, j+1 : j+cols)
    // clang-format off
    kernel_matrix_.block(cur_index * num_params, 
                         cur_index * num_params, 
                         num_params, 
                         num_params) 
                         += weight * ref_kernel;
    // clang-format on
  }
  return true;
}

uint32_t Spline1dKernel::FindIndex(const double x) const {
  auto upper_bound = std::upper_bound(x_knots_.begin() + 1, x_knots_.end(), x);
  return std::min(static_cast<uint32_t>(x_knots_.size() - 1),
                  static_cast<uint32_t>(upper_bound - x_knots_.begin())) -
         1;
}

void Spline1dKernel::AddDistanceOffset(const double weight) {
  for (uint32_t i = 1; i < x_knots_.size(); ++i) {
    const double cur_x = x_knots_[i] - x_knots_[i - 1];
    double       pw_x  = 2.0 * weight;
    for (uint32_t j = 0; j < spline_order_ + 1; ++j) {
      offset_((i - 1) * (spline_order_ + 1) + j, 0) -= pw_x;
      pw_x *= cur_x;
    }
  }
}
}  // namespace planning
}  // namespace apollo
