#include <iostream>
#include <vector>

#include "Eigen/Core"

using namespace std;
Eigen::MatrixXd kernel_fx_;
Eigen::MatrixXd kernel_derivative_;
Eigen::MatrixXd kernel_second_order_derivative_;
Eigen::MatrixXd kernel_third_order_derivative_;

void CalculateFx(const uint32_t num_params) {
  cout << "kernel_fx:" << endl;

  kernel_fx_ = Eigen::MatrixXd::Zero(num_params, num_params);
  for (int r = 0; r < kernel_fx_.rows(); ++r) {
    for (int c = 0; c < kernel_fx_.cols(); ++c) {
      kernel_fx_(r, c) = 1.0 / (r + c + 1.0);
      cout << to_string(1) << "/" << to_string(int(r + c + 1.0)) << "\t";
    }
    cout << endl;
  }
}

void CalculateDerivative(const uint32_t num_params) {
  cout << "kernel_derivative:" << endl;

  kernel_derivative_ = Eigen::MatrixXd::Zero(num_params, num_params);
  for (int r = 0; r < kernel_derivative_.rows(); ++r) {
    for (int c = 0; c < kernel_derivative_.cols(); ++c) {
      if (r < 1 || c < 1) {
        cout << "0\t";
        continue;
      }

      kernel_derivative_(r, c) = r * c / (r + c - 1.0);
      cout << to_string(r * c) << "/" << to_string(int(r + c - 1.0)) << "\t";
    }
    cout << endl;
  }
}

void CalculateSecondOrderDerivative(const uint32_t num_params) {
  cout << "kernel_second_order_derivative:" << endl;

  kernel_second_order_derivative_ = Eigen::MatrixXd::Zero(num_params, num_params);
  for (int r = 0; r < kernel_second_order_derivative_.rows(); ++r) {
    for (int c = 0; c < kernel_second_order_derivative_.cols(); ++c) {
      if (r < 2 || c < 2) {
        cout << "0\t";
        continue;
      }
      kernel_second_order_derivative_(r, c) = (r * r - r) * (c * c - c) / (r + c - 3.0);
      cout << to_string((r * r - r) * (c * c - c)) << "/" << to_string(int(r + c - 3.0)) << "\t";
    }
    cout << endl;
  }
}

void CalculateThirdOrderDerivative(const uint32_t num_params) {
  cout << "kernel_third_order_derivative:" << endl;

  kernel_third_order_derivative_ = Eigen::MatrixXd::Zero(num_params, num_params);
  for (int r = 0; r < kernel_third_order_derivative_.rows(); ++r) {
    for (int c = 0; c < kernel_third_order_derivative_.cols(); ++c) {
      if (r < 3 || c < 3) {
        cout << "0\t";
        continue;
      }

      kernel_third_order_derivative_(r, c) =
          (r * r - r) * (r - 2) * (c * c - c) * (c - 2) / (r + c - 5.0);
      cout << kernel_third_order_derivative_(r, c) << "\t";
    }
    cout << endl;
  }
}

void eigen_block() {
  int                 num_params = 6;
  std::vector<double> power_x    = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  Eigen::MatrixXd ref_kernel(num_params, num_params);
  for (uint32_t r = 0; r < num_params; ++r) {
    for (uint32_t c = 0; c < num_params; ++c) {
      // 可以注意到每个元素前都乘以了2，这是为了和二次优化问题的一般形式中的1/2进行抵消的
      ref_kernel(r, c) = power_x[r + c];
    }
  }
  /*
   *
   *
   * 0       1       2       3       4       5
   * 1       2       3       4       5       6
   * 2       3       4       5       6       7
   * 3       4       5       6       7       8
   * 4       5       6       7       8       9
   * 5       6       7       8       9       10
   *
   */

  for (uint32_t r = 0; r < num_params; ++r) {
    for (uint32_t c = 0; c < num_params; ++c) {
      cout << ref_kernel(r, c) << "\t";
    }
    cout << endl;
  }

  auto block = ref_kernel.block(2, 2, 2, 2);
  //   auto block = ref_kernel.block<3, 3>(1, 1);
  for (uint32_t r = 0; r < 2; ++r) {
    for (uint32_t c = 0; c < 2; ++c) {
      cout << block(r, c) << "\t";
    }
    cout << endl;
  }
}

bool func(bool* arr, int len) {
  bool ret = false;
  for (int i = 0; i < len; i++) {
    ret |= !arr[i];
  }
  return ret;
}

void func11(int x) {
  // int                 x          = 2;
  int                 num_params = 6;
  std::vector<double> x_pow(2 * num_params + 1, 1.0);
  x_pow.reserve(2 * num_params + 1);
  for (uint32_t i = 1; i < 2 * num_params + 1; ++i) {
    x_pow[i] = x_pow[i - 1] * x;
    cout << " " << x_pow[i];
  }
  cout << endl;

  for (uint32_t r = 1; r < num_params; ++r) {
    for (uint32_t c = 1; c < num_params; ++c) {
      cout << " " << x_pow[r + c - 1];
    }
    cout << endl;
  }
}

int main() {
  // CalculateFx(6);
  // CalculateDerivative(6);
  // CalculateSecondOrderDerivative(6);
  // CalculateThirdOrderDerivative(6);
  func11(1);
  func11(2);
  // func11(3);
  return 0;
}