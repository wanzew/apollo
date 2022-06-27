## 螺旋曲线的构造

apollo中的螺旋曲线是弧长-切线角度五次多项式函数：

![[公式]](https://www.zhihu.com/equation?tex=%5Ctheta%3Das%5E%7B5%7D%2Bbs%5E%7B4%7D%2Bcs%5E%7B3%7D%2Bds%5E%7B2%7D%2Bes%2Bf)

由于这种曲线的曲率是弧长的多项式函数，可以在路径规划中很好地控制曲率变化，螺旋曲线的曲率函数为：

![[公式]](https://www.zhihu.com/equation?tex=%5Ckappa%3D%5Cfrac%7Bd%5Ctheta%7D%7Bds%7D%3D5as%5E%7B4%7D%2B4bs%5E%7B3%7D%2B3cs%5E%7B2%7D%2B2ds%2Be)

螺旋曲线类QuinticSpiralPath位置在modules/planning/math/curve1d/quintic_spiral_path.h，其继承了五阶多项式曲线，我们在[Planning 基础库——多项式曲线类](https://zhuanlan.zhihu.com/p/445115282)中讲了仅需要曲线起点0点和终点p的函数值，二阶导函数值和三阶导函数值即可求得系数： ![[公式]](https://www.zhihu.com/equation?tex=f%28x%29%3Dax%5E%7B5%7D%2Bbx%5E%7B4%7D%2Bcx%5E%7B3%7D%2Bdx%5E%7B2%7D%2Bex%2Bf%5C%5C+%5Cleft%5C%7B++++++++++++++++%5Cbegin%7Barray%7D%7B%2A%2Alr%2A%2A%7D++++++++++++++++a%3D%5Cfrac%7B-6f_%7B0%7D%7D%7Bp%5E%7B5%7D%7D-%5Cfrac%7B3f%27_%7B0%7D%7D%7Bp%5E%7B4%7D%7D-%5Cfrac%7B0.5f%27%27_%7B0%7D%7D%7Bp%5E%7B3%7D%7D%2B%5Cfrac%7B6f_%7Bp%7D%7D%7Bp%5E%7B5%7D%7D-%5Cfrac%7B3f%27_%7Bp%7D%7D%7Bp%5E%7B4%7D%7D%2B%5Cfrac%7B0.5f%27%27_%7Bp%7D%7D%7Bp%5E%7B3%7D%7D++%5C%5C++++++++++++++++b%3D%5Cfrac%7B15f_%7B0%7D%7D%7Bp%5E%7B4%7D%7D%2B%5Cfrac%7B8f%27_%7B0%7D%7D%7Bp%5E%7B3%7D%7D%2B%5Cfrac%7B1.5f%27%27_%7B0%7D%7D%7Bp%5E%7B2%7D%7D-%5Cfrac%7B15f_%7Bp%7D%7D%7Bp%5E%7B4%7D%7D%2B%5Cfrac%7B7f%27_%7Bp%7D%7D%7Bp%5E%7B3%7D%7D-%5Cfrac%7Bf%27%27_%7Bp%7D%7D%7Bp%5E%7B2%7D%7D+%5C%5C+++++++++++++++c%3D-%5Cfrac%7B10f_%7B0%7D%7D%7Bp%5E%7B3%7D%7D-%5Cfrac%7B6f%27_%7B0%7D%7D%7Bp%5E%7B2%7D%7D-%5Cfrac%7B1.5f%27%27_%7B0%7D%7D%7Bp%7D%2B%5Cfrac%7B10f_%7Bp%7D%7D%7Bp%5E%7B3%7D%7D-%5Cfrac%7B4f%27_%7Bp%7D%7D%7Bp%5E%7B2%7D%7D%2B%5Cfrac%7B0.5f%27%27_%7Bp%7D%7D%7Bp%7D+%5C%5C+d%3D0.5f%27%27_%7B0%7D+%5C%5C+e%3Df%27_%7B0%7D%5C%5C+f%3Df_%7B0%7D++++++++++++++%5Cend%7Barray%7D+++%5Cright.+++)

```cpp
QuinticSpiralPath(const std::array<double, 3>& start,
                    const std::array<double, 3>& end, const double delta_s);

  QuinticSpiralPath(const double theta0, const double kappa0,
                    const double dkappa0, const double theta1,
                    const double kappa1, const double dkappa1,
                    const double delta_s);
```

## 螺旋曲线笛卡尔坐标变换

螺旋曲线转化为笛卡尔坐标系，x坐标可以将每一段弧长ds投影到x方向ds*cos(θ)，y坐标将每一段弧长ds投影到y方向ds*sin(θ)，加上曲线初始状态的x，y值积分得到

![[公式]](https://www.zhihu.com/equation?tex=x%28s%29%3Dx_%7B0%7D%2B%5Cint_%7B0%7D%5E%7Bs%7Dcos%28%5Ctheta%28s%29%29ds%5C%5C+y%28s%29%3Dy_%7B0%7D%2B%5Cint_%7B0%7D%5E%7Bs%7Dsin%28%5Ctheta%28s%29%29ds)

我们构造螺旋曲线通常都会让初始值x0,y0设置为0，代入高斯积分得到

![[公式]](https://www.zhihu.com/equation?tex=x%28s%29%3D%5Cfrac%7Bs%7D%7B2%7D%5Csum_%7Bi%3D1%7D%5E%7Bn%7D%7Bcos%28%5Ctheta%28%5Cfrac%7Bs%7D%7B2%7D%5Cxi_i%2B%5Cfrac%7Bs%7D%7B2%7D%29%29%7D%5C%5C++y%28s%29%3D%5Cfrac%7Bs%7D%7B2%7D%5Csum_%7Bi%3D1%7D%5E%7Bn%7D%7Bsin%28%5Ctheta%28%5Cfrac%7Bs%7D%7B2%7D%5Cxi_i%2B%5Cfrac%7Bs%7D%7B2%7D%29%29%7D)

从函数中可以看出，如果需要计算x，y坐标值需要对非线性函数进行积分操作，在[Planning 基础库——积分方法](https://zhuanlan.zhihu.com/p/445031799)中我们讲到了几种积分方法，这里使用了高斯-勒让德进行积分求解

对应代码函数接口:

```cpp
template <size_t N>
  double ComputeCartesianDeviationX(const double s) const {
    auto cos_theta = [this](const double s) {
      const auto a = Evaluate(0, s);
      return std::cos(a);
    };
    return common::math::IntegrateByGaussLegendre<N>(cos_theta, 0.0, s);
  }

  template <size_t N>
  double ComputeCartesianDeviationY(const double s) const {
    auto sin_theta = [this](const double s) {
      const auto a = Evaluate(0, s);
      return std::sin(a);
    };
    return common::math::IntegrateByGaussLegendre<N>(sin_theta, 0.0, s);
  }
```

## 偏导数矩阵

我们在观察QuinticSpiralPath类的定义中可以发现，除了具有五次多项式类QuinticPolynomialCurve1d的coef_多项式系数外，还多了一个系数矩阵coef_deriv_

```cpp
  std::array<std::array<double, 7>, 6> coef_deriv_;
```

这个系数矩阵是干嘛的呢？进一步分析我们看到，用到这个系数矩阵的主要是几个尾缀为Derivative的函数

```cpp
 std::pair<double, double> DeriveCartesianDeviation(const size_t param_index)

  double DeriveKappaDerivative(const size_t param_index,
                               const double ratio) const;

  double DeriveDKappaDerivative(const size_t param_index,
                                const double ratio) const;

  double DeriveD2KappaDerivative(const size_t param_index,
                                 const double r) const;
```

这几个函数显然是求导数的，但是我们求导数有了多项式类中Evaluate()了吗，可以求各阶导数的函数值，为什么螺旋曲线类还有求导数的函数呢

```cpp
double Evaluate(const std::uint32_t order, const double p) const override;
```

其实这几个函数是用于在螺旋曲线平滑过程中求偏导数用的，详细可以参见()介绍

螺旋曲线平滑算法中定义了5个优化变量，

![img](https://pic1.zhimg.com/80/v2-4988eb9a79851e840428c2c549680fc4_720w.png)

对于每一条螺旋曲线，需要分别对θ(0),θ'(0),θ''(0),θ(Δs),θ'(Δs),θ'(Δs)求偏导数。刚才我们介绍了五次多项式可以通过初始点状态和终止点构造出来。所以我们尝试对θ(0)求偏导

![[公式]](https://www.zhihu.com/equation?tex=%5Cfrac%7B%5Cdelta+f%28s%29%7D%7Bd%5Ctheta_%7B0%7D%7D%3D++%3D+%5Cfrac+%7B-6s%5E%7B5%7D%7D%7Bp%5E%7B5%7D%7D%2B%5Cfrac+%7B15s%5E%7B4%7D%7D%7Bp%5E%7B4%7D%7D%2B%5Cfrac+%7B-10s%5E%7B3%7D%7D%7Bp%5E%7B3%7D%7D%2B1%5C%5C+%5Cfrac%7B%5Cdelta+f%28s%29%7D%7Bds%27_%7B0%7D%7D%3D%5Cfrac+%7B-3s%5E%7B5%7D%7D%7Bp%5E%7B4%7D%7D%2B%5Cfrac+%7B8s%5E%7B4%7D%7D%7Bp%5E%7B3%7D%7D%2B%5Cfrac+%7B-6s%5E%7B3%7D%7D%7Bp%5E%7B2%7D%7D%2B1%5C%5C+%5C%5C...)

所以我们可以看出来，coef_deriv_系数矩阵正好对应着a,b,c,d,e,f各次参数的矩阵去掉θ(0),θ'(0),θ''(0),θ(Δs),θ'(Δs),θ''(Δs)的参数值，也就是可以归纳为下式（这里遵从apollo中矩阵顺序从低次到高次排序，令p=Δs）

![[公式]](https://www.zhihu.com/equation?tex=f%28%5Ctheta%29%3DA%5Ccdot+%5Cbegin%7Bpmatrix%7D+1%5C%5C++%5Ctheta%5C%5C++%5Ctheta%5E%7B2%7D%5C%5C++%5Ctheta%5E%7B3%7D%5C%5C++%5Ctheta%5E%7B4%7D%5C%5C++%5Ctheta%5E%7B5%7D%5C%5C+%5Cend%7Bpmatrix%7D%5C%5C+A%3D%5Cbegin%7Bpmatrix%7D+%5Ctheta_%7B0%7D+%26+%5Ctheta%27_%7B0%7D+%26+%5Ctheta%27%27_%7B0%7D+%26+%5Ctheta_%7Bp%7D+%26+%5Ctheta%27_%7Bp%7D+%26+%5Ctheta%27%27_%7Bp%7D+%5Cend%7Bpmatrix%7DC_%7Bderiv%7D%5C%5C++C_%7Bderiv%7D%3D%5Cbegin%7BBmatrix%7D++++1%26++0%26++0%26++0%26+0+%260+%5C%5C+++0%26++1%26++0%26+0+%26+0+%260+%5C%5C+++0%26++0%26++0.5%26++0%260++%260+%5C%5C+++%5Cfrac%7B-10%7D%7Bp%5E%7B3%7D%7D%26++%5Cfrac%7B-6%7D%7Bp%5E%7B2%7D%7D%26+%5Cfrac%7B-1.5%7D%7Bp%7D%26+%5Cfrac%7B10%7D%7Bp%5E%7B3%7D%7D%26%5Cfrac%7B-4%7D%7Bp%5E%7B2%7D%7D%26%5Cfrac%7B0.5%7D%7Bp%7D+%5C%5C+++%5Cfrac%7B15%7D%7Bp%5E%7B4%7D%7D%26+%5Cfrac%7B8%7D%7Bp%5E%7B3%7D%7D+%26++%5Cfrac%7B1.5%7D%7Bp%5E%7B2%7D%7D%26%5Cfrac%7B-15%7D%7Bp%5E%7B4%7D%7D++%26+%5Cfrac%7B7%7D%7Bp%5E%7B3%7D%7D+%26%5Cfrac%7B-1%7D%7Bp%5E%7B2%7D%7D+%5C%5C+%5Cfrac%7B-6%7D%7Bp%5E%7B5%7D%7D%26+%5Cfrac%7B-3%7D%7Bp%5E%7B4%7D%7D+%26+%5Cfrac%7B-0.5%7D%7Bp%5E%7B3%7D%7D++%26+%5Cfrac%7B6%7D%7Bp%5E%7B5%7D%7D++%26+%5Cfrac%7B-3%7D%7Bp%5E%7B4%7D%7D++%26+%5Cfrac%7B-0.5%7D%7Bp%5E%7B3%7D%7D+%5Cend%7BBmatrix%7D)

代码中coef_deriv_矩阵是7x6的矩阵，当前矩阵只是6x6的矩阵，那么第7列是什么呢。

我们的优化变量中还有一个Δs，也就是式中的p，优化过程中我们还需要将系数a,b,c,d,e,f对p求偏导

![[公式]](https://www.zhihu.com/equation?tex=C_%7Bp%7D+%3D+%5Cbegin%7Bpmatrix%7D+%5Cfrac%7B%5Cdelta+f%7D%7Bdp%7D%5C%5C+%5Cfrac%7B%5Cdelta+e%7D%7Bdp%7D%5C%5C+%5Cfrac%7B%5Cdelta+d%7D%7Bdp%7D%5C%5C+%5Cfrac%7B%5Cdelta+c%7D%7Bdp%7D%5C%5C+%5Cfrac%7B%5Cdelta+b%7D%7Bdp%7D%5C%5C+%5Cfrac%7B%5Cdelta+a%7D%7Bdp%7D%5C%5C+%5Cend%7Bpmatrix%7D+%3D+%5Cbegin%7Bpmatrix%7D+0%5C%5C++0%5C%5C++0%5C%5C++%5Cfrac%7B30%5Ctheta_%7B0%7D%7D%7Bp%5E%7B4%7D%7D%2B%5Cfrac%7B12%5Ctheta%27_%7B0%7D%7D%7Bp%5E%7B3%7D%7D%2B%5Cfrac%7B1.5%5Ctheta%27%27_%7B0%7D%7D%7Bp%5E%7B2%7D%7D%2B%5Cfrac%7B30%5Ctheta_%7Bp%7D%7D%7Bp%5E%7B4%7D%7D%2B%5Cfrac%7B8%5Ctheta%27_%7Bp%7D%7D%7Bp%5E%7B3%7D%7D-%5Cfrac%7B0.5%5Ctheta%27%27_%7Bp%7D%7D%7Bp%5E%7B2%7D%7D%5C%5C++%5Cfrac%7B-60%5Ctheta_%7B0%7D%7D%7Bp%5E%7B5%7D%7D-%5Cfrac%7B24%5Ctheta%27_%7B0%7D%7D%7Bp%5E%7B4%7D%7D-%5Cfrac%7B3%5Ctheta%27%27_%7B0%7D%7D%7Bp%5E%7B3%7D%7D%2B%5Cfrac%7B60%5Ctheta_%7Bp%7D%7D%7Bp%5E%7B5%7D%7D-%5Cfrac%7B21%5Ctheta%27_%7Bp%7D%7D%7Bp%5E%7B4%7D%7D%2B%5Cfrac%7B2%5Ctheta%27%27_%7Bp%7D%7D%7Bp%5E%7B3%7D%7D%5C%5C++%5Cfrac%7B30%5Ctheta_%7B0%7D%7D%7Bp%5E%7B6%7D%7D%2B%5Cfrac%7B12%5Ctheta%27_%7B0%7D%7D%7Bp%5E%7B5%7D%7D%2B%5Cfrac%7B1.5%5Ctheta%27%27_%7B0%7D%7D%7Bp%5E%7B4%7D%7D-%5Cfrac%7B30%5Ctheta_%7Bp%7D%7D%7Bp%5E%7B6%7D%7D%2B%5Cfrac%7B12%5Ctheta%27_%7Bp%7D%7D%7Bp%5E%7B5%7D%7D-%5Cfrac%7B1.5%5Ctheta%27%27_%7Bp%7D%7D%7Bp%5E%7B4%7D%7D+%5Cend%7Bpmatrix%7D)

## 偏导函数推导

假如我们对于曲线上的点si求对于p的偏导 ![[公式]](https://www.zhihu.com/equation?tex=s_%7Bi%7D%3Dp%5Ccdot+r_i%5C%5C+%5Cfrac%7B%5Cdelta+f%28s_%7Bi%7D%29%7D%7Bdp%7D%3D+%5Cfrac%7B%5Cdelta+a%7D%7Bdp%7D+s_%7Bi%7D%5E%7B5%7D%2B5s_%7Bi%7D%5E4r_%7Bi%7D%2B+%5Cfrac%7B%5Cdelta+b%7D%7Bdp%7D+s_%7Bi%7D%5E%7B4%7D%2B4s_%7Bi%7D%5E3r_%7Bi%7D%2B+%5Cfrac%7B%5Cdelta+c%7D%7Bdp%7D+s_%7Bi%7D%5E%7B3%7D%2B3s_%7Bi%7D%5E2r_%7Bi%7D%2B+%5Cfrac%7B%5Cdelta+d%7D%7Bdp%7D+s_%7Bi%7D%5E%7B2%7D%2B2s_%7Bi%7Dr_%7Bi%7D%2B+%5Cfrac%7B%5Cdelta+e%7D%7Bdp%7D+s_%7Bi%7D%2Br_%7Bi%7D%2B+%5Cfrac%7B%5Cdelta+f%7D%7Bdp%7D+)

对应函数接口如下

```cpp
double QuinticSpiralPath::DeriveTheta(const size_t param_index,
 const double r) const {
 double s = param_ * r;
 double s2 = s * s;
 double s3 = s2 * s;
 double s4 = s2 * s2;
 double s5 = s3 * s2;

 double derivative =
 coef_deriv_[5][param_index] * s5 + coef_deriv_[4][param_index] * s4 +
 coef_deriv_[3][param_index] * s3 + coef_deriv_[2][param_index] * s2 +
 coef_deriv_[1][param_index] * s + coef_deriv_[0][param_index];

 if (param_index == DELTA_S) {
 derivative += coef_[5] * 5.0 * s4 * r + coef_[4] * 4.0 * s3 * r +
 coef_[3] * 3.0 * s2 * r + coef_[2] * 2.0 * s * r +
 coef_[1] * r;
  }
 return derivative;
}
```

其中param_index表示求偏导的变量，如下所示，0和1分别表示了螺旋曲线的起点和终点

```cpp
  static const size_t THETA0 = 0;
  static const size_t KAPPA0 = 1;
  static const size_t DKAPPA0 = 2;
  static const size_t THETA1 = 3;
  static const size_t KAPPA1 = 4;
  static const size_t DKAPPA1 = 5;
  static const size_t DELTA_S = 6;
```

优化中我们还需要对螺旋曲线的曲率函数求偏导

曲率函数我们可以得到

![[公式]](https://www.zhihu.com/equation?tex=%5Ckappa%28s%29%3D+%5Cfrac%7B%5Ctheta%28s%29%7D%7Bds%7D+%3D5as%5E4%2B4bs%5E3%2B3cs%5E2%2B2ds%2Be)

对a,b,c,d,e求导，带入C_deriv矩阵可以得到θ(0),θ'(0),θ''(0),θ(Δs),θ'(Δs),θ'(Δs)的偏导数

![[公式]](https://www.zhihu.com/equation?tex=%5Cbegin%7Bpmatrix%7D+%5Cfrac%7B%5Ckappa+%28s%29%7D%7Bd%5Ctheta_0%7D%5C%5C+%5Cfrac%7B%5Ckappa+%28s%29%7D%7Bd%5Ckappa_0%7D%5C%5C+%5Cfrac%7B%5Ckappa+%28s%29%7D%7Bdd%5Ckappa_0%7D%5C%5C+%5Cfrac%7B%5Ckappa+%28s%29%7D%7Bd%5Ctheta_1%7D%5C%5C+%5Cfrac%7B%5Ckappa+%28s%29%7D%7Bd%5Ckappa_1%7D%5C%5C+%5Cfrac%7B%5Ckappa+%28s%29%7D%7Bdd%5Ckappa_1%7D%5C%5C+%5Cend%7Bpmatrix%7D%3D+%5Cbegin%7Bpmatrix%7D+0+%261++%262++%26+3+%264++%265++%5Cend%7Bpmatrix%7D+C_%7Bderiv%7D+%5Cbegin%7Bpmatrix%7D+1%5C%5C+s%5C%5C+s%5E2%5C%5C+s%5E3%5C%5C+s%5E4%5C%5C+s%5E5++%5Cend%7Bpmatrix%7D)

对应代码接口

```cpp
double QuinticSpiralPath::DeriveKappaDerivative(const size_t param_index,
                                                const double r) const {
  double s = param_ * r;
  double s2 = s * s;
  double s3 = s2 * s;
  double s4 = s2 * s2;

  double derivative = 5.0 * coef_deriv_[5][param_index] * s4 +
                      4.0 * coef_deriv_[4][param_index] * s3 +
                      3.0 * coef_deriv_[3][param_index] * s2 +
                      2.0 * coef_deriv_[2][param_index] * s +
                      coef_deriv_[1][param_index];

  if (param_index == DELTA_S) {
    derivative += 5.0 * coef_[5] * 4.0 * s3 * r +
                  4.0 * coef_[4] * 3.0 * s2 * r + 3.0 * coef_[3] * 2.0 * s * r +
                  2.0 * coef_[2] * r;
  }
  return derivative;
}
```

同理，曲率变化率的偏导函数可以表示为：

![[公式]](https://www.zhihu.com/equation?tex=%5Ckappa%28s%29%3D+%5Cfrac%7B%5Ctheta%28s%29%7D%7Bds%7D+%3D20as%5E3%2B12bs%5E2%2B6cs%5E1%2B2d)

偏导函数：

![[公式]](https://www.zhihu.com/equation?tex=%5Cbegin%7Bpmatrix%7D+%5Cfrac%7B%5Ckappa+%28s%29%7D%7Bd%5Ctheta_0%7D%5C%5C+%5Cfrac%7B%5Ckappa+%28s%29%7D%7Bd%5Ckappa_0%7D%5C%5C+%5Cfrac%7B%5Ckappa+%28s%29%7D%7Bdd%5Ckappa_0%7D%5C%5C+%5Cfrac%7B%5Ckappa+%28s%29%7D%7Bd%5Ctheta_1%7D%5C%5C+%5Cfrac%7B%5Ckappa+%28s%29%7D%7Bd%5Ckappa_1%7D%5C%5C+%5Cfrac%7B%5Ckappa+%28s%29%7D%7Bdd%5Ckappa_1%7D%5C%5C+%5Cend%7Bpmatrix%7D%3D+%5Cbegin%7Bpmatrix%7D+0+%260++%262++%26+6+%2612++%2620++%5Cend%7Bpmatrix%7D+C_%7Bderiv%7D+%5Cbegin%7Bpmatrix%7D+1%5C%5C+s%5C%5C+s%5E2%5C%5C+s%5E3%5C%5C+s%5E4%5C%5C+s%5E5++%5Cend%7Bpmatrix%7D)

对应代码接口

```cpp
double QuinticSpiralPath::DeriveDKappaDerivative(const size_t param_index,
                                                 const double r) const {
  double s = param_ * r;
  double s2 = s * s;
  double s3 = s2 * s;

  double derivative = 20.0 * coef_deriv_[5][param_index] * s3 +
                      12.0 * coef_deriv_[4][param_index] * s2 +
                      6.0 * coef_deriv_[3][param_index] * s +
                      2.0 * coef_deriv_[2][param_index];

  if (param_index == DELTA_S) {
    derivative += 20.0 * coef_[5] * 3.0 * s2 * r +
                  12.0 * coef_[4] * 2.0 * s * r + 6.0 * coef_[3] * r;
  }
  return derivative;
}
```

## 笛卡尔坐标的偏导数

我们还需要计算x，y相对于θ(0),θ'(0),θ''(0),θ(Δs),θ'(Δs),θ'(Δs)的偏导数：

已知

![[公式]](https://www.zhihu.com/equation?tex=x%28s%29%3D%5Cfrac%7Bs%7D%7B2%7D%5Csum_%7Bi%3D1%7D%5E%7Bn%7D%7Bcos%28%5Ctheta%28%5Cfrac%7Bs%7D%7B2%7D%5Cxi_i%2B%5Cfrac%7Bs%7D%7B2%7D%29%29%7D%5C%5C++y%28s%29%3D%5Cfrac%7Bs%7D%7B2%7D%5Csum_%7Bi%3D1%7D%5E%7Bn%7D%7Bsin%28%5Ctheta%28%5Cfrac%7Bs%7D%7B2%7D%5Cxi_i%2B%5Cfrac%7Bs%7D%7B2%7D%29%29%7D)

分别对x(s),y(s)求对于θ(0)的偏导

![[公式]](https://www.zhihu.com/equation?tex=%5Cfrac%7Bd+x%7D%7Bd%5Ctheta_0%7D%3D+%5Cfrac%7Bd+x%7D%7Bd%5Ctheta%7D+%5Cfrac%7Bd+%5Ctheta%7D%7Bd%5Ctheta_0%7D%3D+-%5Cfrac%7Bs%7D%7B2%7D+%5Csum_%7Bi%3D1%7D%5E%7Bn%7D%7Bsin%28%5Ctheta%28%5Cfrac%7Bs%7D%7B2%7D%5Cxi_i%2B%5Cfrac%7Bs%7D%7B2%7D%29%29%7D%5Ccdot+%5Cfrac%7Bd+%5Ctheta%28%5Cfrac%7Bs%7D%7B2%7D%5Cxi_i%2B%5Cfrac%7Bs%7D%7B2%7D%29%7D%7Bd%5Ctheta_0%7D%5C%5C++%5Cfrac%7Bd+y%7D%7Bd%5Ctheta_0%7D%3D+%5Cfrac%7Bd+y%7D%7Bd%5Ctheta%7D+%5Cfrac%7Bd+%5Ctheta%7D%7Bd%5Ctheta_0%7D%3D+%5Cfrac%7Bs%7D%7B2%7D+%5Csum_%7Bi%3D1%7D%5E%7Bn%7D%7Bcos%28%5Ctheta%28%5Cfrac%7Bs%7D%7B2%7D%5Cxi_i%2B%5Cfrac%7Bs%7D%7B2%7D%29%29%7D%5Ccdot+%5Cfrac%7Bd+%5Ctheta%28%5Cfrac%7Bs%7D%7B2%7D%5Cxi_i%2B%5Cfrac%7Bs%7D%7B2%7D%29%7D%7Bd%5Ctheta_0%7D%5C%5C+)

其中θ(s)对θ(0)的偏导数前文已经求出。同理可以推导出x,y对θ'(0),θ''(0),θ(Δs),θ'(Δs),θ'(Δs)的偏导数

对应代码接口如下：

```cpp
template <size_t N>
  std::pair<double, double> DeriveCartesianDeviation(
      const size_t param_index) const {
    auto gauss_points = common::math::GetGaussLegendrePoints<N>();
    std::array<double, N> x = gauss_points.first;
    std::array<double, N> w = gauss_points.second;

    std::pair<double, double> cartesian_deviation = {0.0, 0.0};
    for (size_t i = 0; i < N; ++i) {
      double r = 0.5 * x[i] + 0.5;
      auto curr_theta = Evaluate(0, r * param_);
      double derived_theta = DeriveTheta(param_index, r);

      cartesian_deviation.first +=
          w[i] * (-std::sin(curr_theta)) * derived_theta;
      cartesian_deviation.second += w[i] * std::cos(curr_theta) * derived_theta;
    }

    cartesian_deviation.first *= param_ * 0.5;
    cartesian_deviation.second *= param_ * 0.5;

    if (param_index == DELTA_S) {
      for (size_t i = 0; i < N; ++i) {
        double r = 0.5 * x[i] + 0.5;
        auto theta_angle = Evaluate(0, r * param_);

        cartesian_deviation.first += 0.5 * w[i] * std::cos(theta_angle);
        cartesian_deviation.second += 0.5 * w[i] * std::sin(theta_angle);
      }
    }
    return cartesian_deviation;
  }
```