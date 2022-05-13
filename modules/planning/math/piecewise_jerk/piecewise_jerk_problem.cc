/******************************************************************************
 * Copyright 2018 The Apollo Authors. All Rights Reserved.
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

#include "modules/planning/math/piecewise_jerk/piecewise_jerk_problem.h"

#include "cyber/common/log.h"
#include "modules/planning/common/planning_gflags.h"

namespace apollo {
namespace planning {

namespace {
constexpr double kMaxVariableRange = 1.0e10;
}  // namespace

PiecewiseJerkProblem::PiecewiseJerkProblem(const size_t                 num_of_knots,
                                           const double                 delta_s,
                                           const std::array<double, 3>& x_init) {
  CHECK_GE(num_of_knots, 2U);
  num_of_knots_ = num_of_knots;
  x_init_       = x_init;
  delta_s_      = delta_s;
  x_bounds_.resize(num_of_knots_, std::make_pair(-kMaxVariableRange, kMaxVariableRange));
  dx_bounds_.resize(num_of_knots_, std::make_pair(-kMaxVariableRange, kMaxVariableRange));
  ddx_bounds_.resize(num_of_knots_, std::make_pair(-kMaxVariableRange, kMaxVariableRange));
  weight_x_ref_vec_ = std::vector<double>(num_of_knots_, 0.0);
}

// FormulateProblem 这个函数用于构造最优化问题具体矩阵。首先构造出P矩阵即代价函数，
// 然后构造A矩阵即约束矩阵以及上下边界lower_bounds和upper_bounds，最后构建一次项q向量。
// 构造完后将矩阵都存储进OSQPData这个结构体里，以便后续直接调用osqp库进行求解。
// 以下详细解读这些矩阵的具体形式。
OSQPData* PiecewiseJerkProblem::FormulateProblem() {
  // calculate kernel
  std::vector<c_float> P_data;
  std::vector<c_int>   P_indices;
  std::vector<c_int>   P_indptr;

  // CalculateKnernel
  // 用于构造代价函数即P矩阵，由于代码过长就不全部贴出来了。我在初次看这块代码的时候，
  // 就被csc_matrix这奇葩的矩阵构造方式折磨了好久，并且我看网上许多讲解Apollo
  // 二次规划的文章也都没有具体到矩阵的实际形式的。我这里把代价函数P矩阵实际解构了出来，
  // 我觉得对于第一次接触这块代码的朋友，能够直观的看到这个矩阵还是会对理解Apollo的算法思想有很大帮助的。
  // 可以注意到每个元素前都乘以了2，这是为了和二次优化问题的一般形式中的1/2进行抵消的。
  CalculateKernel(&P_data, &P_indices, &P_indptr);

  // calculate affine constraints
  std::vector<c_float> A_data;
  std::vector<c_int>   A_indices;
  std::vector<c_int>   A_indptr;
  std::vector<c_float> lower_bounds;
  std::vector<c_float> upper_bounds;
  CalculateAffineConstraint(&A_data, &A_indices, &A_indptr, &lower_bounds, &upper_bounds);

  // calculate offset
  std::vector<c_float> q;
  CalculateOffset(&q);

  OSQPData* data = reinterpret_cast<OSQPData*>(c_malloc(sizeof(OSQPData)));
  CHECK_EQ(lower_bounds.size(), upper_bounds.size());

  size_t kernel_dim            = 3 * num_of_knots_;
  size_t num_affine_constraint = lower_bounds.size();

  data->n = kernel_dim;
  data->m = num_affine_constraint;
  data->P = csc_matrix(kernel_dim, kernel_dim, P_data.size(), CopyData(P_data), CopyData(P_indices),
                       CopyData(P_indptr));
  data->q = CopyData(q);
  data->A = csc_matrix(num_affine_constraint, kernel_dim, A_data.size(), CopyData(A_data),
                       CopyData(A_indices), CopyData(A_indptr));
  data->l = CopyData(lower_bounds);
  data->u = CopyData(upper_bounds);
  return data;
}

bool PiecewiseJerkProblem::Optimize(const int max_iter) {
  // 该函数中会调用FormulateProblem来构造出二次规划问题的框架，再调用osqp库进行求解，从而求出最优path
  // 需要注意的是，二次规划问题的求解方式有许多种，包括拉格朗日法，梯度下降等等，Apollo
  // 采用的osqp这个第三方库个人感觉实际用时求解效率还是比较高的，也比较好用。
  // 唯一不太舒服的地方在于其矩阵的构造形式为csc_matrix，该种方法构造矩阵不够直观，比较复杂，
  // 后面我会把矩阵用LaTex打出来，让读者可以更直观理解。
  OSQPData* data = FormulateProblem();

  OSQPSettings* settings = SolverDefaultSettings();
  settings->max_iter     = max_iter;

  OSQPWorkspace* osqp_work = nullptr;
  // osqp_work                = osqp_setup(data, settings);
  osqp_setup(&osqp_work, data, settings);

  osqp_solve(osqp_work);

  auto status = osqp_work->info->status_val;

  if (status < 0 || (status != 1 && status != 2)) {
    AERROR << "failed optimization status:\t" << osqp_work->info->status;
    osqp_cleanup(osqp_work);
    FreeData(data);
    c_free(settings);
    return false;
  } else if (osqp_work->solution == nullptr) {
    AERROR << "The solution from OSQP is nullptr";
    osqp_cleanup(osqp_work);
    FreeData(data);
    c_free(settings);
    return false;
  }

  // extract primal results
  x_.resize(num_of_knots_);
  dx_.resize(num_of_knots_);
  ddx_.resize(num_of_knots_);
  for (size_t i = 0; i < num_of_knots_; ++i) {
    x_.at(i)   = osqp_work->solution->x[i] / scale_factor_[0];
    dx_.at(i)  = osqp_work->solution->x[i + num_of_knots_] / scale_factor_[1];
    ddx_.at(i) = osqp_work->solution->x[i + 2 * num_of_knots_] / scale_factor_[2];
  }

  // Cleanup
  osqp_cleanup(osqp_work);
  FreeData(data);
  c_free(settings);
  return true;
}

// 该函数是用于构造最优化问题的限制条件，即A矩阵的。需要注意的是，该方法为基类中实现，因此对于path优化和speed优化来说，
// 两者在Apollo 的算法思想里，所收到的约束是一样的。Apollo
// 该算法的精妙之处就在于，将path和speed分别在SL和ST空间中进行考虑，使得两者的优化思想非常类似，
// 很巧妙地完成两个维度的求解。但与此同时，我感觉这也限制了speed优化，对于动态障碍物的处理就不够完备，后面我会单独再写文章详细讲这一块。
void PiecewiseJerkProblem::CalculateAffineConstraint(std::vector<c_float>* A_data,
                                                     std::vector<c_int>*   A_indices,
                                                     std::vector<c_int>*   A_indptr,
                                                     std::vector<c_float>* lower_bounds,
                                                     std::vector<c_float>* upper_bounds) {
  // 3N params bounds on x, x', x''
  // 3(N-1) constraints on x, x', x''
  // 3 constraints on x_init_
  const int n                  = static_cast<int>(num_of_knots_);
  const int num_of_variables   = 3 * n;
  const int num_of_constraints = num_of_variables + 3 * (n - 1) + 3;
  lower_bounds->resize(num_of_constraints);
  upper_bounds->resize(num_of_constraints);

  std::vector<std::vector<std::pair<c_int, c_float>>> variables(num_of_variables);

  int constraint_index = 0;
  // set x, x', x'' bounds
  for (int i = 0; i < num_of_variables; ++i) {
    if (i < n) {
      variables[i].emplace_back(constraint_index, 1.0);
      lower_bounds->at(constraint_index) = x_bounds_[i].first * scale_factor_[0];
      upper_bounds->at(constraint_index) = x_bounds_[i].second * scale_factor_[0];
    } else if (i < 2 * n) {
      variables[i].emplace_back(constraint_index, 1.0);

      lower_bounds->at(constraint_index) = dx_bounds_[i - n].first * scale_factor_[1];
      upper_bounds->at(constraint_index) = dx_bounds_[i - n].second * scale_factor_[1];
    } else {
      variables[i].emplace_back(constraint_index, 1.0);
      lower_bounds->at(constraint_index) = ddx_bounds_[i - 2 * n].first * scale_factor_[2];
      upper_bounds->at(constraint_index) = ddx_bounds_[i - 2 * n].second * scale_factor_[2];
    }
    ++constraint_index;
  }
  CHECK_EQ(constraint_index, num_of_variables);

  // x(i->i+1)''' = (x(i+1)'' - x(i)'') / delta_s
  for (int i = 0; i + 1 < n; ++i) {
    variables[2 * n + i].emplace_back(constraint_index, -1.0);
    variables[2 * n + i + 1].emplace_back(constraint_index, 1.0);
    lower_bounds->at(constraint_index) = dddx_bound_.first * delta_s_ * scale_factor_[2];
    upper_bounds->at(constraint_index) = dddx_bound_.second * delta_s_ * scale_factor_[2];
    ++constraint_index;
  }

  // x(i+1)' - x(i)' - 0.5 * delta_s * x(i)'' - 0.5 * delta_s * x(i+1)'' = 0
  for (int i = 0; i + 1 < n; ++i) {
    variables[n + i].emplace_back(constraint_index, -1.0 * scale_factor_[2]);
    variables[n + i + 1].emplace_back(constraint_index, 1.0 * scale_factor_[2]);
    variables[2 * n + i].emplace_back(constraint_index, -0.5 * delta_s_ * scale_factor_[1]);
    variables[2 * n + i + 1].emplace_back(constraint_index, -0.5 * delta_s_ * scale_factor_[1]);
    lower_bounds->at(constraint_index) = 0.0;
    upper_bounds->at(constraint_index) = 0.0;
    ++constraint_index;
  }

  // x(i+1) - x(i) - delta_s * x(i)'
  // - 1/3 * delta_s^2 * x(i)'' - 1/6 * delta_s^2 * x(i+1)''
  auto delta_s_sq_ = delta_s_ * delta_s_;
  for (int i = 0; i + 1 < n; ++i) {
    variables[i].emplace_back(constraint_index, -1.0 * scale_factor_[1] * scale_factor_[2]);
    variables[i + 1].emplace_back(constraint_index, 1.0 * scale_factor_[1] * scale_factor_[2]);
    variables[n + i].emplace_back(constraint_index,
                                  -delta_s_ * scale_factor_[0] * scale_factor_[2]);
    variables[2 * n + i].emplace_back(constraint_index,
                                      -delta_s_sq_ / 3.0 * scale_factor_[0] * scale_factor_[1]);
    variables[2 * n + i + 1].emplace_back(constraint_index,
                                          -delta_s_sq_ / 6.0 * scale_factor_[0] * scale_factor_[1]);

    lower_bounds->at(constraint_index) = 0.0;
    upper_bounds->at(constraint_index) = 0.0;
    ++constraint_index;
  }

  // constrain on x_init
  variables[0].emplace_back(constraint_index, 1.0);
  lower_bounds->at(constraint_index) = x_init_[0] * scale_factor_[0];
  upper_bounds->at(constraint_index) = x_init_[0] * scale_factor_[0];
  ++constraint_index;

  variables[n].emplace_back(constraint_index, 1.0);
  lower_bounds->at(constraint_index) = x_init_[1] * scale_factor_[1];
  upper_bounds->at(constraint_index) = x_init_[1] * scale_factor_[1];
  ++constraint_index;

  variables[2 * n].emplace_back(constraint_index, 1.0);
  lower_bounds->at(constraint_index) = x_init_[2] * scale_factor_[2];
  upper_bounds->at(constraint_index) = x_init_[2] * scale_factor_[2];
  ++constraint_index;

  CHECK_EQ(constraint_index, num_of_constraints);

  int ind_p = 0;
  for (int i = 0; i < num_of_variables; ++i) {
    A_indptr->push_back(ind_p);
    for (const auto& variable_nz : variables[i]) {
      // coefficient
      A_data->push_back(variable_nz.second);

      // constraint index
      A_indices->push_back(variable_nz.first);
      ++ind_p;
    }
  }
  // We indeed need this line because of
  // https://github.com/oxfordcontrol/osqp/blob/master/src/cs.c#L255
  A_indptr->push_back(ind_p);
}

OSQPSettings* PiecewiseJerkProblem::SolverDefaultSettings() {
  // Define Solver default settings
  OSQPSettings* settings = reinterpret_cast<OSQPSettings*>(c_malloc(sizeof(OSQPSettings)));
  osqp_set_default_settings(settings);
  settings->polish             = true;
  settings->verbose            = FLAGS_enable_osqp_debug;
  settings->scaled_termination = true;
  return settings;
}

void PiecewiseJerkProblem::set_x_bounds(std::vector<std::pair<double, double>> x_bounds) {
  CHECK_EQ(x_bounds.size(), num_of_knots_);
  x_bounds_ = std::move(x_bounds);
}

void PiecewiseJerkProblem::set_dx_bounds(std::vector<std::pair<double, double>> dx_bounds) {
  CHECK_EQ(dx_bounds.size(), num_of_knots_);
  dx_bounds_ = std::move(dx_bounds);
}

void PiecewiseJerkProblem::set_ddx_bounds(std::vector<std::pair<double, double>> ddx_bounds) {
  CHECK_EQ(ddx_bounds.size(), num_of_knots_);
  ddx_bounds_ = std::move(ddx_bounds);
}

void PiecewiseJerkProblem::set_x_bounds(const double x_lower_bound, const double x_upper_bound) {
  for (auto& x : x_bounds_) {
    x.first  = x_lower_bound;
    x.second = x_upper_bound;
  }
}

void PiecewiseJerkProblem::set_dx_bounds(const double dx_lower_bound, const double dx_upper_bound) {
  for (auto& x : dx_bounds_) {
    x.first  = dx_lower_bound;
    x.second = dx_upper_bound;
  }
}

void PiecewiseJerkProblem::set_ddx_bounds(const double ddx_lower_bound,
                                          const double ddx_upper_bound) {
  for (auto& x : ddx_bounds_) {
    x.first  = ddx_lower_bound;
    x.second = ddx_upper_bound;
  }
}

void PiecewiseJerkProblem::set_x_ref(const double weight_x_ref, std::vector<double> x_ref) {
  CHECK_EQ(x_ref.size(), num_of_knots_);
  weight_x_ref_ = weight_x_ref;
  // set uniform weighting
  weight_x_ref_vec_ = std::vector<double>(num_of_knots_, weight_x_ref);
  x_ref_            = std::move(x_ref);
  has_x_ref_        = true;
}

void PiecewiseJerkProblem::set_x_ref(std::vector<double> weight_x_ref_vec,
                                     std::vector<double> x_ref) {
  CHECK_EQ(x_ref.size(), num_of_knots_);
  CHECK_EQ(weight_x_ref_vec.size(), num_of_knots_);
  // set piecewise weighting
  weight_x_ref_vec_ = std::move(weight_x_ref_vec);
  x_ref_            = std::move(x_ref);
  has_x_ref_        = true;
}

void PiecewiseJerkProblem::set_end_state_ref(const std::array<double, 3>& weight_end_state,
                                             const std::array<double, 3>& end_state_ref) {
  weight_end_state_  = weight_end_state;
  end_state_ref_     = end_state_ref;
  has_end_state_ref_ = true;
}

void PiecewiseJerkProblem::FreeData(OSQPData* data) {
  delete[] data->q;
  delete[] data->l;
  delete[] data->u;

  delete[] data->P->i;
  delete[] data->P->p;
  delete[] data->P->x;

  delete[] data->A->i;
  delete[] data->A->p;
  delete[] data->A->x;
}

}  // namespace planning
}  // namespace apollo
