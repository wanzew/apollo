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
 * @file : spline_1d_generator.cc
 * @brief: piecewise_smoothing_spline (pss) generator class
 *           solve pss by qp algorithm, include adding constraint, adding
 *kernel, and solver solve
 **/

#include "modules/planning/math/smoothing_spline/spline_1d_generator.h"

#include <algorithm>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "modules/common/log.h"
#include "modules/common/math/qp_solver/active_set_qp_solver.h"
#include "modules/common/math/qp_solver/qp_solver_gflags.h"
#include "modules/common/time/time.h"
#include "modules/planning/common/planning_gflags.h"

namespace apollo {
namespace planning {
namespace {

constexpr double kMaxBound = 1e3;
}

using apollo::common::time::Clock;
using Eigen::MatrixXd;

Spline1dGenerator::Spline1dGenerator(const std::vector<double>& x_knots,
                                     const uint32_t             spline_order)
    : spline_(x_knots, spline_order)
    , spline_constraint_(x_knots, spline_order)
    , spline_kernel_(x_knots, spline_order) {}

void Spline1dGenerator::Reset(const std::vector<double>& x_knots, const uint32_t spline_order) {
  spline_            = Spline1d(x_knots, spline_order);
  spline_constraint_ = Spline1dConstraint(x_knots, spline_order);
  spline_kernel_     = Spline1dKernel(x_knots, spline_order);
}

Spline1dConstraint* Spline1dGenerator::mutable_spline_constraint() { return &spline_constraint_; }

Spline1dKernel* Spline1dGenerator::mutable_spline_kernel() { return &spline_kernel_; }

// 求解：
// 由于qpoase所需的输入都是一维的，因此需要对矩阵、各约束进行转换。以下就对求解过程中的各个参数的赋值和转换进行了分析。

bool Spline1dGenerator::Solve() {
  // 二维hissen矩阵，定义为MatrixXd::Zero(total_params_, total_params_)
  // 线性矩阵ｇ，实则一个向量　定义为MatrixXd::Zero(total_params_, 1);
  auto ineq_constraint_matrix   = spline_constraint_.inequality_constraint().constraint_matrix();
  auto ineq_constraint_boundary = spline_constraint_.inequality_constraint().constraint_boundary();
  auto eq_constraint_matrix     = spline_constraint_.equality_constraint().constraint_matrix();
  auto eq_constraint_boundary   = spline_constraint_.equality_constraint().constraint_boundary();

  const MatrixXd& kernel_matrix                  = spline_kernel_.kernel_matrix();
  const MatrixXd& offset                         = spline_kernel_.offset();
  const MatrixXd& inequality_constraint_matrix   = ineq_constraint_matrix;
  const MatrixXd& inequality_constraint_boundary = ineq_constraint_boundary;
  const MatrixXd& equality_constraint_matrix     = eq_constraint_matrix;
  const MatrixXd& equality_constraint_boundary   = eq_constraint_boundary;

  if (kernel_matrix.rows() != kernel_matrix.cols()) {
    AERROR << "kernel_matrix.rows() [" << kernel_matrix.rows() << "] and kernel_matrix.cols() ["
           << kernel_matrix.cols() << "] should be identical.";
    return false;
  }

  // clang-format off
  int num_param      = kernel_matrix.rows();
  int num_constraint = equality_constraint_matrix.rows()  
                     + inequality_constraint_matrix.rows();
  bool use_hotstart  = last_problem_success_ &&
                      (FLAGS_enable_sqp_solver && 
                       sqp_solver_ != nullptr &&
                       num_param == last_num_param_ && 
                       num_constraint == last_num_constraint_);
  //clang-format on

  // use_hot_start：当上一次有解且和这一次的约束 参数数量相同就不再进行初始化和init设置了
  if (!use_hotstart) {
    // 1. qpoase参数设置
    //设置变量个数 约束个数 Hessian矩阵类型未知
    sqp_solver_.reset(new ::qpOASES::SQProblem(num_param, num_constraint, ::qpOASES::HST_UNKNOWN));
    ::qpOASES::Options my_options;
    //设定Specifies the frequency of full refactorisation of proj. Hessian
    my_options.enableCholeskyRefactorisation = 1;

    my_options.epsNum = FLAGS_default_active_set_eps_num;  // 设定 比率测试的分子公差。
    my_options.epsDen = FLAGS_default_active_set_eps_den;  // 设定比率测试的分母公差。
    my_options.epsIterRef = FLAGS_default_active_set_eps_iter_ref;  // 设定迭代优化的早期终止公差
    sqp_solver_->setOptions(my_options);                            // 设定参数选项
    if (!FLAGS_default_enable_active_set_debug_info) {              // 设定不输出
      sqp_solver_->setPrintLevel(qpOASES::PL_NONE);
    }
  }

  // definition of qpOASESproblem
  const int kNumOfMatrixElements = kernel_matrix.rows() * kernel_matrix.cols();
  double    h_matrix[kNumOfMatrixElements];

  const int kNumOfOffsetRows = offset.rows();
  double    g_matrix[kNumOfOffsetRows];  // ， 与offset一般大
  int       index = 0;

  // kernel_matrix.rows() == kernel_matrix.cols() == offset.rows()
  for (int r = 0; r < kernel_matrix.rows(); ++r) {
    g_matrix[r] = offset(r, 0);
    for (int c = 0; c < kernel_matrix.cols(); ++c) {
      h_matrix[index++] = kernel_matrix(r, c);  // h矩阵为逐行赋值
    }
  }
  DCHECK_EQ(index, kernel_matrix.rows() * kernel_matrix.cols());

  // search space lower bound and uppper bound
  double lower_bound[num_param];
  double upper_bound[num_param];

  const double l_lower_bound_ = -kMaxBound;
  const double l_upper_bound_ = kMaxBound;
  for (int i = 0; i < num_param; ++i) {
    lower_bound[i] = l_lower_bound_;
    upper_bound[i] = l_upper_bound_;
  }

  // constraint matrix construction
  double affine_constraint_matrix[num_param * num_constraint];  // 大小为参数个数乘以约束个数
  double constraint_lower_bound[num_constraint];                // 大小为约束条件个数
  double constraint_upper_bound[num_constraint];                // 大小为约束条件个数

  index = 0;
  for (int r = 0; r < equality_constraint_matrix.rows(); ++r) {
    constraint_lower_bound[r] = equality_constraint_boundary(r, 0);  // 上下限相同
    constraint_upper_bound[r] = equality_constraint_boundary(r, 0);

    for (int c = 0; c < num_param; ++c) {
      // 按照行顺序输入矩阵A中
      affine_constraint_matrix[index++] = equality_constraint_matrix(r, c);
    }
  }

  DCHECK_EQ(index, equality_constraint_matrix.rows() * num_param);

  // 增加不等式约束
  const double constraint_upper_bound_ = kMaxBound;
  for (int r = 0; r < inequality_constraint_matrix.rows(); ++r) {
    constraint_lower_bound[r + equality_constraint_boundary.rows()] =  // 下限由不等式获取
        inequality_constraint_boundary(r, 0);
    // 上限是最大值30，已经将所有上边界进行转化为下边界了
    constraint_upper_bound[r + equality_constraint_boundary.rows()] =  //
        constraint_upper_bound_;

    for (int c = 0; c < num_param; ++c) {
      // 按照行顺序输入矩阵A中
      affine_constraint_matrix[index++] = inequality_constraint_matrix(r, c);
    }
  }
  DCHECK_EQ(index, equality_constraint_matrix.rows() * num_param +
                inequality_constraint_boundary.rows() * num_param);

  // initialize problem
  int max_iteration_ = 1000;
  int max_iter       = std::max(max_iteration_, num_constraint);

  ::qpOASES::returnValue ret;
  const double           start_timestamp = Clock::NowInSeconds();
  if (use_hotstart) {
    ADEBUG << "using SQP hotstart.";
    // 初始化和求解
    // clang-format off
    // 如果能够热启动，使用部分上一次的搜索结果
    ret = sqp_solver_->hotstart(h_matrix,
                                g_matrix, 
                                affine_constraint_matrix, 
                                lower_bound,
                                upper_bound, 
                                constraint_lower_bound, 
                                constraint_upper_bound,
                                max_iter);
    if (ret != qpOASES::SUCCESSFUL_RETURN) {
      AERROR << "Fail to hotstart spline 1d, will use re-init instead.";
      // 重新求解
      ret = sqp_solver_->init(h_matrix, 
                              g_matrix, 
                              affine_constraint_matrix, 
                              lower_bound, 
                              upper_bound,
                              constraint_lower_bound, 
                              constraint_upper_bound, 
                              max_iter);
    }
  } else {
    ADEBUG << "no using SQP hotstart.";
    ret = sqp_solver_->init(h_matrix, 
                            g_matrix, 
                            affine_constraint_matrix, 
                            lower_bound, 
                            upper_bound,
                            constraint_lower_bound, 
                            constraint_upper_bound, 
                            max_iter);
    // clang-format on
  }
  const double end_timestamp = Clock::NowInSeconds();
  ADEBUG << "Spline1dGenerator QP solve time: " << (end_timestamp - start_timestamp) * 1000
         << " ms.";

  if (ret != qpOASES::SUCCESSFUL_RETURN) {
    if (ret == qpOASES::RET_MAX_NWSR_REACHED) {
      AERROR << "qpOASES solver failed due to reached max iteration";
    } else {
      AERROR << "qpOASES solver failed due to infeasibility or other internal "
                "reasons:"
             << ret;
    }
    last_problem_success_ = false;
    return false;
  }

  last_problem_success_ = true;
  // 3. 获取结果
  double result[num_param];                //
  memset(result, 0, sizeof result);        // 全为0
  sqp_solver_->getPrimalSolution(result);  // 获取结果

  MatrixXd solved_params = MatrixXd::Zero(num_param, 1);
  for (int i = 0; i < num_param; ++i) {  // 结果转化为内部需要的一维格式
    solved_params(i, 0) = result[i];
    ADEBUG << "spline 1d solved param[" << i << "]: " << result[i];
  }

  last_num_param_      = num_param;
  last_num_constraint_ = num_constraint;

  // 求取参数化的样条曲线
  return spline_.SetSplineSegs(solved_params, spline_.spline_order());
}

const Spline1d& Spline1dGenerator::spline() const { return spline_; }

}  // namespace planning
}  // namespace apollo
