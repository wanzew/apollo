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
 * @file : affine_constraint.cc
 **/

#include "modules/planning/math/smoothing_spline/affine_constraint.h"

#include "modules/common/log.h"

namespace apollo {
namespace planning {

AffineConstraint::AffineConstraint(const bool is_equality)
    : is_equality_(is_equality) {}

AffineConstraint::AffineConstraint(const Eigen::MatrixXd& constraint_matrix,
                                   const Eigen::MatrixXd& constraint_boundary,
                                   const bool             is_equality)
    : constraint_matrix_(constraint_matrix)
    , constraint_boundary_(constraint_boundary)
    , is_equality_(is_equality) {
  CHECK_EQ(constraint_boundary.rows(), constraint_matrix.rows());
}

void AffineConstraint::SetIsEquality(const double is_equality) { is_equality_ = is_equality; }

const Eigen::MatrixXd& AffineConstraint::constraint_matrix() const { return constraint_matrix_; }

const Eigen::MatrixXd& AffineConstraint::constraint_boundary() const {
  return constraint_boundary_;
}

//追加约束，是最底层的函数
bool AffineConstraint::AddConstraint(const Eigen::MatrixXd& constraint_matrix,
                                     const Eigen::MatrixXd& constraint_boundary) {
  // 假设添加的约束是 Ax=b， constraint_matrix 就是A，constraint_boundary 就是b，对5次多项式就是
  //[1,p,p^2,p^3,p^4,p^5]，p是多项式曲线上某点的横坐标
  // x就是未知的多项式系数
  // 自然A的一行对应的Ax的计算结果就是b的一行，所以 A.rows==b.rows。
  // b的元素就是目标方程的计算结果，是一个实数，自然 b.cols==1
  if (constraint_matrix.rows() != constraint_boundary.rows()) {
    AERROR << "Fail to add constraint because constraint matrix rows != "
              "constraint boundary rows.";
    AERROR << "constraint matrix rows = " << constraint_matrix.rows();
    AERROR << "constraint boundary rows = " << constraint_boundary.rows();
    return false;
  }

  if (constraint_matrix_.rows() == 0) {
    constraint_matrix_   = constraint_matrix;
    constraint_boundary_ = constraint_boundary;
    return true;
  }
  if (constraint_matrix_.cols() != constraint_matrix.cols()) {
    AERROR << "constraint_matrix_ cols and constraint_matrix cols do not match.";
    AERROR << "constraint_matrix_.cols() = " << constraint_matrix_.cols();
    AERROR << "constraint_matrix.cols() = " << constraint_matrix.cols();
    return false;
  }
  if (constraint_boundary.cols() != 1) {
    AERROR << "constraint_boundary.cols() should be 1.";
    return false;
  }

  Eigen::MatrixXd n_matrix(constraint_matrix_.rows() + constraint_matrix.rows(),
                           constraint_matrix_.cols());
  Eigen::MatrixXd n_boundary(constraint_boundary_.rows() + constraint_boundary.rows(), 1);

  // 追加约束
  n_matrix << constraint_matrix_, constraint_matrix;
  n_boundary << constraint_boundary_, constraint_boundary;
  constraint_matrix_   = n_matrix;
  constraint_boundary_ = n_boundary;
  return true;
}

}  // namespace planning
}  // namespace apollo
