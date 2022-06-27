# 分段五次螺旋线非线性优化拟合(基于Apollo)

[TOC]
分段五次螺旋线非线性优化拟合(基于Apollo)



## 螺旋线介绍

曲线方程可以表示为以弧长为参数的参数表达式, 好处是在使用弧长参数表达式下的 CnCn 连续性转换成任何参数形式仍然能保证 CnCn 连续, 同时如果表达为 θ(s)θ(s) 即切线的弧长参数表达式可以很方便地求得曲率的 nn 阶导.

## 问题模型建立

如下图
[![illustraion.jpg](https://www.chuxin911.com/piecewise_cubic_spiral_curve_20210923/illustraion.jpg)](https://www.chuxin911.com/piecewise_cubic_spiral_curve_20210923/illustraion.jpg)

[illustraion.jpg](https://www.chuxin911.com/piecewise_cubic_spiral_curve_20210923/illustraion.jpg)


给定 nn 个离散点 P0(x0,y0)P0(x0,y0)~Pn−1(xn−1,yn−1)Pn−1(xn−1,yn−1), 被分为 n−1n−1 段, 每段长度为 ΔsiΔsi.
每段曲线的参数表达式如下:



θi(s)=ais5+bis4+cis3+dis2+eis+fi→θi(s)=→pi⋅→sθi(s)=ais5+bis4+cis3+dis2+eis+fi→θi(s)=pi→⋅s→

其中 s∈[0,Δsi]s∈[0,Δsi].

解析求得两点确定的一条曲线(第 ii 段)过程如下:
未知参数向量 →pipi→ 为 6 行向量: ai,bi,ci,di,ei,fiai,bi,ci,di,ei,fi, 需要 6 个方程求解得出.

θi(0)=θiθi(Δsi)=θi+1˙θi(0)=˙θi˙θi(Δsi)=˙θi+1¨θi(0)=¨θiθi(Δsi)=¨θi+1θi(0)=θiθi(Δsi)=θi+1θi˙(0)=θi˙θi˙(Δsi)=θi+1˙θi¨(0)=θi¨θi(Δsi)=θi+1¨

其中 θi(0),˙θi(0),¨θi(0)θi(0),θi˙(0),θi¨(0) 分别为曲线起点处的方向, 曲率, 曲率导数, 对于 i+1i+1 时的三个参数分别为曲线终点的方向, 曲率, 曲率导数.
可以使用逆矩阵求解, 例如 Eigen 里的 `reverse()` 等函数. 可以提前解出用参数表达式的话, 提高计算效率, 如下:

ai=−6θis5−3˙θis4−¨θi2s3+6θi+1s5−3˙θi+1s4+¨θi+12s3bi=15θis4+8˙θis3+3¨θi2s2−15θi+1s4+7˙θi+1s3−¨θi+1s2ci=−10θis3−6˙θis2−3¨θi2s+10θi+1s3−4˙θi+1s2+¨θi+12sdi=¨θi2;ei=˙θifi=θiai=−6θis5−3θi˙s4−θi¨2s3+6θi+1s5−3θi+1˙s4+θi+1¨2s3bi=15θis4+8θi˙s3+3θi¨2s2−15θi+1s4+7θi+1˙s3−θi+1¨s2ci=−10θis3−6θi˙s2−3θi¨2s+10θi+1s3−4θi+1˙s2+θi+1¨2sdi=θi¨2;ei=θi˙fi=θi

为了保证相邻段曲线之间的连续性 (C2)(C2) 连续, 可得如下方程:

xi+1=xi+∫Δsi0cos(θ(s))dsyi+1=yi+∫Δsi0sin(θ(s))dsθi+1=θi(Δsi)˙θi+1=˙θi(Δsi)¨θi+1=¨θi(Δsi)xi+1=xi+∫0Δsicos⁡(θ(s))dsyi+1=yi+∫0Δsisin⁡(θ(s))dsθi+1=θi(Δsi)θi+1˙=θi˙(Δsi)θi+1¨=θi¨(Δsi)

## 优化模型

### 优化变量

优化变量设置为 →qq→ 拟合所需要的结果, 这样不需要再用系数求一遍, 同时如果使用系数作为优化变量的话, 无法直接在约束中添加曲率/曲率变化率的bounds, 需要构造 constraints 或者是间接添加 bounds. 如下每个元素均为 nn 行的向量.

→q=[→θ,→˙θ,→¨θ,→x,→y,→Δs]q→=[θ→,θ˙→,θ¨→,x→,y→,Δs→]例如→Δs=[Δs0Δs1⋯Δsn−1]Δs→=[Δs0Δs1⋯Δsn−1]其中Δs0=0Δs0=0

### 目标方程

fcost=wlengthn−1∑i=0Δsi+wkppan−1∑i=0(m−1∑j=0(˙θij(s))2)+wdkppan−1∑i=0(m−1∑j=0(¨θij(s))2)fcost=wlength∑i=0n−1Δsi+wkppa∑i=0n−1(∑j=0m−1(θij˙(s))2)+wdkppa∑i=0n−1(∑j=0m−1(θij¨(s))2)

分别代表长度, 曲率, 曲率变化率, mm 表示对第 ii 段曲线长度等分为 mm 子段(分段的原因, 初步猜测为减小数值积分误差).

### 约束条件

包含 bounds 与 constraints, 分别用 BB 与 CC 序号表示, 下表带 ii 表示共 nn 个.

#### 起点等式约束(bounds)

θ0=θstart⋯(B0)κ0=κstart⋯(B1)˙κ0=˙κstart⋯(B2)x0=xref0⋯(B3)y0=yref0⋯(B4)θ0=θstart⋯(B0)κ0=κstart⋯(B1)κ0˙=κstart˙⋯(B2)x0=xref0⋯(B3)y0=yref0⋯(B4)

#### 终点等式约束(bounds)

θn−1=θend⋯(B5)κn−1=κend⋯(B6)˙κn−1=˙κend⋯(B7)xn−1=xrefn−1⋯(B8)yn−1=yrefn−1⋯(B9)θn−1=θend⋯(B5)κn−1=κend⋯(B6)κn−1˙=κend˙⋯(B7)xn−1=xrefn−1⋯(B8)yn−1=yrefn−1⋯(B9)

#### 物理模型约束(bounds)

最大转角/转角变化率/转角变化率的变化率( Apollo 中标定值).

θi−1−π2≤θi≤θi−1+π2⋯(B10i)−0.25≤κ≤+0.25⋯(B11i)−0.02≤˙κ≤+0.02⋯(B12i)θi−1−π2≤θi≤θi−1+π2⋯(B10i)−0.25≤κ≤+0.25⋯(B11i)−0.02≤κ˙≤+0.02⋯(B12i)

#### 中间点范围约束(bounds)

规定连续两个点之间最大转向为四分之一圆弧.

xrefi−ri≤xi≤xrefi+ri⋯(B13i)yrefi−ri≤yi≤yrefi+ri⋯(B14i)Di−2ri≤Δsi≤Diπ2⋯(B15i)其中,Di=√(xrefi−xrefi+1)2+(yrefi−yrefi+1)2xrefi−ri≤xi≤xrefi+ri⋯(B13i)yrefi−ri≤yi≤yrefi+ri⋯(B14i)Di−2ri≤Δsi≤Diπ2⋯(B15i)其中,Di=(xrefi−xrefi+1)2+(yrefi−yrefi+1)2

#### 连接点等式约束(bounds+constraints)

xi+1=xi+∫Δsi0cos(θ(s))ds⋯(C0i)yi+1=yi+∫Δsi0sin(θ(s))ds⋯(C1i)θi+1=θi(Δsi)⋯(B16i)˙θi+1=˙θi(Δsi)⋯(B17i)¨θi+1=¨θi(Δsi)⋯(B18i)xi+1=xi+∫0Δsicos⁡(θ(s))ds⋯(C0i)yi+1=yi+∫0Δsisin⁡(θ(s))ds⋯(C1i)θi+1=θi(Δsi)⋯(B16i)θi+1˙=θi˙(Δsi)⋯(B17i)θi+1¨=θi¨(Δsi)⋯(B18i)

#### 位置平移非等式约束(constraints)

(xi−xrefi)2+(yi−yrefi)2≤r2i⋯(C2i)(xi−xrefi)2+(yi−yrefi)2≤ri2⋯(C2i)

### 数值积分求解弧长

连接点等式约束中出现了积分, 其中存在 Fresnel 积分, 求不出解析解, Apollo 中采用 Gauss-Legendre 积分方法. 采用此积分方法的好处是计算效率较高, 只需计算特定几个点处的函数值加权和即可逼近积分值, 并且能保证较满意的精度. Apollo 中采用的是五点高斯积分.
下面为 Gauss-Legendre 积分公式
∫1−1f(x)dx≈n∑k=0Akf(xk)∫−11f(x)dx≈∑k=0nAkf(xk)
下图为 Gauss-Legendre 积分积分点与权重系数表:
[![Gauss_Legendre.png](https://www.chuxin911.com/piecewise_cubic_spiral_curve_20210923/Gauss_Legendre.png)](https://www.chuxin911.com/piecewise_cubic_spiral_curve_20210923/Gauss_Legendre.png)

[Gauss_Legendre.png](https://www.chuxin911.com/piecewise_cubic_spiral_curve_20210923/Gauss_Legendre.png)



具体步骤:

1. 需要将自变量区间 [0,Δsi][0,Δsi] 线性变换为 t∗[−1,1]t∗[−1,1].
2. 在 n=4n=4 处的 xixi 点处求函数值的权重 AiAi 和.
   具体过程参考 Apollo 代码: `modules/common/math/interal.cc`(math 里面有很多小工具值得借鉴一下.)

### 求解库

Apollo 采用 IPOPT 非线性求解器求解此问题, IPOPT 的使用指南可以参考[博文](http://www.chuxin911.com/IPOPT_intro_20210906/). 只需要将接口所需的参数/表达式传递进去即可. 下面的部分为具体的操作, Apollo 的代码会在后面逐行解析头文件.

#### 起点: Ipopt::TNLP::get_starting_point

把 θ0=θstart,κ0=κstart,˙κ0=˙κstart,x0=xref0,y0=yref0θ0=θstart,κ0=κstart,κ0˙=κstart˙,x0=xref0,y0=yref0 传进去.

#### 目标方程:Ipopt::TNLP::eval_f

把 fcostfcost 传进去.

#### 目标方程Gradient向量: Ipopt::TNLP::eval_grad_f

这个是整个求解过程中较为复杂的一部分, 此部分求偏导的结果可以用于约束方程的偏导求解. 优化变量为 6n6n 向量, 因此目标方程的 Gradient 向量也为 6n6n 向量, 拿任意相邻 2 组 12 行为例, [θi,˙θi,¨θi,xi,yi,Δsi,θi+1,˙θi+1,¨θi+1,xi+1,yi+1,Δsi+1]T[θi,θi˙,θi¨,xi,yi,Δsi,θi+1,θi+1˙,θi+1¨,xi+1,yi+1,Δsi+1]T. fcostfcost 中只需考虑求和公式中的第 i,i+1i,i+1 项即可, 其余项对优化变量的偏导均为 0. 注, 后面有时会省略下标 ii, 偷懒一下 :-). 同时对 ΔsΔs 与 ss 的求导结果一致.
即:

fcost(i,i+1)=wlength[Δsi(→q)+Δsi+1(→q)]+wkppa[(˙θi(→q))2+(˙θi+1(→q))2]+wdkppa[(¨θi(→q))2+(¨θi+1(→q))2]fcost(i,i+1)=wlength[Δsi(q→)+Δsi+1(q→)]+wkppa[(θi˙(q→))2+(θi+1˙(q→))2]+wdkppa[(θi¨(q→))2+(θi+1¨(q→))2]

- 对 →x,→yx→,y→ 的偏导为 0.∂˙θi2/∂˙θi+12/∂¨θi2/∂¨θi+12/∂Δsi∂xi/∂yi/∂xi+1/∂yi+1=∂fcost(i,i+1)∂xi/∂yi/∂xi+1/∂yi+1=0∂θi˙2/∂θi+1˙2/∂θi¨2/∂θi+1¨2/∂Δsi∂xi/∂yi/∂xi+1/∂yi+1=∂fcost(i,i+1)∂xi/∂yi/∂xi+1/∂yi+1=0
- 对 θi,˙θi,¨θi,θi+1,˙θi+1,¨θi+1θi,θi˙,θi¨,θi+1,θi+1˙,θi+1¨ 求偏导可以参考如下例子:∂˙θi(→q)2∂θi=m∑r=02˙θi(s)∗∂(ais5+bis4+cis3+dis2+eis+fi)∂θi∣∣∣s=rj=m∑r=02˙θi(s)∗∂[(−6θiΔs5−3˙θiΔs4−¨θi2Δs3+6θi+1Δs5−3˙θi+1Δs4+¨θi+12Δs3)s5+bis4+cis3+dis2+eis+fi]∂θi∣∣∣s=rj=m∑r=02˙θi(rj)∗(5−6Δs5r4j+415Δs4r3j+3−10Δs3r2j)∂θi(q→)˙2∂θi=∑r=0m2θi(s)˙∗∂(ais5+bis4+cis3+dis2+eis+fi)∂θi|s=rj=∑r=0m2θi(s)˙∗∂[(−6θiΔs5−3θi˙Δs4−θi¨2Δs3+6θi+1Δs5−3θi+1˙Δs4+θi+1¨2Δs3)s5+bis4+cis3+dis2+eis+fi]∂θi|s=rj=∑r=0m2θi(rj)˙∗(5−6Δs5rj4+415Δs4rj3+3−10Δs3rj2)其中rj=ratioj∗Δsratioj=j/mrj=ratioj∗Δsratioj=j/m同理将参数向量 →qiq→i 带入可得 ¨θi(→q)2θi(q→)¨2 的偏导.
  另外∂Δs∂θ/∂˙θ/∂¨θ=0∂Δs∂θ/∂θ˙/∂θ¨=0
- 对 Δsi,Δsi+1Δsi,Δsi+1 求偏导∂Δsi(→q)∂Δsi=1∂˙θi(→q)2∂Δsi=m∑r=02˙θi(s)∗[∂[(−6θiΔs5−3˙θiΔs4−¨θi2Δs3+6θi+1Δs5−3˙θi+1Δs4+¨θi+12Δs3)s5+bis4+cis3+dis2+eis+fi]∂Δs+∂rj∂Δs]∣∣∣s=rj=m∑r=02˙θi(s)∗[(30θiΔs6+12˙θiΔs5+¨3θi2Δs3+−30θi+1Δs6+12˙θi+1Δs5+−3¨θi+12Δs4)s5+∂(bis4+cis3+dis2+eis+fi)+∂rj∂Δs]∣∣∣s=rj∂Δsi(q→)∂Δsi=1∂θi(q→)˙2∂Δsi=∑r=0m2θi(s)˙∗[∂[(−6θiΔs5−3θi˙Δs4−θi¨2Δs3+6θi+1Δs5−3θi+1˙Δs4+θi+1¨2Δs3)s5+bis4+cis3+dis2+eis+fi]∂Δs+∂rj∂Δs]|s=rj=∑r=0m2θi(s)˙∗[(30θiΔs6+12θi˙Δs5+3θi¨2Δs3+−30θi+1Δs6+12θi+1˙Δs5+−3θi+1¨2Δs4)s5+∂(bis4+cis3+dis2+eis+fi)+∂rj∂Δs]|s=rj其中∂rj∂Δs=∂˙θ(ratioj∗s)∂Δs∣∣∣s=rj=ratioj∂˙θ(s)∂s∣∣∣s=rj=ratioj¨θ(s)∣∣∣s=rj∂rj∂Δs=∂θ(ratioj∗s)˙∂Δs|s=rj=ratioj∂θ(s)˙∂s|s=rj=ratiojθ(s)¨|s=rj其余处偏导均为 0. 此处对 ΔsΔs 求偏导出现了加法, 为 Apollo 中计算公式, 笔者认为是不是应该采取乘法.下面的约束方程对 ΔsΔs 求偏导也是如此.

#### 边界: Ipopt::TNLP::get_bounds_info

共 9n+109n+10 个
起终点边界: B0-B9,10 个.
物理模型约束边界: B10-B12,3n 个.
中间点范围边界: B13-B15,3n 个.
连接点等式约束: B16-B18,3n 个.

#### 约束方程: Ipopt::TNLP::eval_g

连接点等式约束:

g1i(→p)=(xi+1−xi−∫Δsi0cos(θ(s))ds)2=0g2i(→p)=(yi+1−yi−∫Δsi0sin(θ(s))ds)2=0g1i(p→)=(xi+1−xi−∫0Δsicos⁡(θ(s))ds)2=0g2i(p→)=(yi+1−yi−∫0Δsisin⁡(θ(s))ds)2=0位置平移非等式约束:g3i(→p)=(xi−xrefi)2+(yi−yrefi)2−r2i≤0g3i(p→)=(xi−xrefi)2+(yi−yrefi)2−ri2≤0

共 3n3n 个.

#### 约束方程: Ipopt::TNLP::eval_jac_g

令 xrj=(xi+1−xi−∫Δsi0cos(θ(s))ds)2∣∣∣s=rjxrj=(xi+1−xi−∫0Δsicos⁡(θ(s))ds)2|s=rj 将积分空间 [0,Δsi][0,Δsi] 转换到 [−1,1][−1,1], 自变量变为 (Δsi2r+Δsi2)(Δsi2r+Δsi2), 以便进行数值积分求解. 对于 g1ig1i:∂g1i∂θi=−2xrjΔsi2(∫1−1sinθ∂θi(→q)∂θi)∣∣∣s=rj∂g1i∂˙θi=−2xrjΔsi2(∫1−1sinθ∂θi(→q)∂˙θi)∣∣∣s=rj∂g1i∂˙θi=−2xrjΔsi2(∫1−1sinθ∂θi(→q)∂¨θi)∣∣∣s=rj∂g1i∂θi+1=−2xrjΔsi2(∫1−1sinθ∂θi(→q)∂θi+1)∣∣∣s=rj∂g1i∂˙θi+1=−2xrjΔsi2(∫1−1sinθ∂θi(→q)∂˙θi+1)∣∣∣s=rj∂g1i∂˙θi+1=−2xrjΔsi2(∫1−1sinθ∂θi(→q)∂¨θi+1)∣∣∣s=rj∂g1i∂xi=−2xrj∗(−1)∂g1i∂xi+1=−2xrj∗(1)∂g1i∂yi=0∂g1i∂yi+1=0∂g1i∂Δsi=∂g1i∂θi+(−2xrj)12∫1−1cos(θ)ds∣∣∣s=rj∂g1i∂Δsi+1=0∂g1i∂θi=−2xrjΔsi2(∫−11sin⁡θ∂θi(q→)∂θi)|s=rj∂g1i∂θi˙=−2xrjΔsi2(∫−11sin⁡θ∂θi(q→)∂θi˙)|s=rj∂g1i∂θi˙=−2xrjΔsi2(∫−11sin⁡θ∂θi(q→)∂θi¨)|s=rj∂g1i∂θi+1=−2xrjΔsi2(∫−11sin⁡θ∂θi(q→)∂θi+1)|s=rj∂g1i∂θi+1˙=−2xrjΔsi2(∫−11sin⁡θ∂θi(q→)∂θi+1˙)|s=rj∂g1i∂θi+1˙=−2xrjΔsi2(∫−11sin⁡θ∂θi(q→)∂θi+1¨)|s=rj∂g1i∂xi=−2xrj∗(−1)∂g1i∂xi+1=−2xrj∗(1)∂g1i∂yi=0∂g1i∂yi+1=0∂g1i∂Δsi=∂g1i∂θi+(−2xrj)12∫−11cos⁡(θ)ds|s=rj∂g1i∂Δsi+1=0相同地可以得到 g2ig2i 的偏导. 对 g3ig3i 仅在 xixi 与 yiyi 处有偏导.∂g3i∂xi=2(xi−xrefi)∂g3i∂yi=2(yi−yrefi)∂g3i∂xi=2(xi−xrefi)∂g3i∂yi=2(yi−yrefi)

#### 目标方程Hessian矩阵: Ipopt::TNLP::eval_h

Hessain 矩阵计算公式如下:

H=σf∇2f(→p)+3n∑i=1λi∇2gi(→p)H=σf∇2f(p→)+∑i=13nλi∇2gi(p→)Apollo 中将 Hessian 矩阵设置为 0, 实际理论上不为零,下面举出 2 个例子: 例1:∂2˙θi(→q)∂θi∂Δsi=m∑r=0∂(5−6Δs5r4j+415Δs4r3j+3−10Δs3r2j)∂Δs=m∑r=0(150Δs6r4j+−240Δs5r3j+90Δs4r2j)∂2θi(q→)˙∂θi∂Δsi=∑r=0m∂(5−6Δs5rj4+415Δs4rj3+3−10Δs3rj2)∂Δs=∑r=0m(150Δs6rj4+−240Δs5rj3+90Δs4rj2)例2:∂2g3i∂2xi=2∂2g3i∂2yi=2∂2g3i∂2xi=2∂2g3i∂2yi=2

设置为 0 可能是考虑到求解性? 由于笔者未跑实际的代码调试, 不敢妄下结论, 后续实际调试 Apollo 的时候, 观察按照理论值求解结果如何.

## Apollo(v6.0)源码解析

### 对所有相关代码逐行解析后, 主要类的 `.h` 文件分析如下

- curve1d 基类, 纯虚类

  ```
  
  ```

- polynomial_curve1d, 继承 1 级类, 纯虚类

  ```
  
  ```

- quintic_polynomial_curve1d, 继承 2 级类, 五次多项式曲线

  ```
  
  ```

- quintic_spiral_path, 继承 3 级类, 五次螺旋线

  ```
  
  ```

- QuinticSpiralPathWithDerivation, 继承 3 级类, 带导数信息的五次螺旋线
  与 quintic_spiral_path 的区别是用 unordered_map 存储了积分点处的 cache_cartesian_deriv_, cache_kappa_deriv_, cache_dkappa_deriv_, 并且仅有头文件.

  ```
  
  ```

- SpiralProblemInterface, 继承Ipopt::TNLP, 使用 IPOPT 的接口

  ```
  
  ```

- PiecewiseQuinticSpiralPath, 继承 1 级类
  把多段定义成 path, 并定义操作, 实际未被使用.

  ```
  
  ```

### 观码感触

- 抽象层次清晰合理:1d曲线->1d多项式曲线->1d五次多项式曲线->五次多项式螺旋.
- 通过继承父类的构造函数构造子类的同时实现多态构造子类构造函数, `QuinticSpiralPathWithDerivation` 类继承 `QuinticPolynomialCurve1d` 类的构造函数.
- 类内构造函数对另一个构造函数的调用.
- 对于计算重复使用的项, 用散列表存储, 提升效率, 例如 `cache_kappa_deriv_`, 散列函数的设计是使用 nn( 最大类型个数)进制.
- 把很多矩阵运算改成了方程式运算, 虽然难懂了一些, 官方称效率能提升.
- 使用 lambda 函数计算 cos 值, 简洁.
- 将多阶的求导/取值/积分利用矩阵存储, 简洁.
- 吐槽: 长度为啥不用 length/l 这种一看就理解的命名, 非要用param_/p/delta_s 乱七八糟的不统一命名.
- 吐槽: theta0, kappa0, dkappa0, theta1, kappa1, dkappa1 的命名刚开始很难懂, 如果命名为 theta_i<->theta_i_1 可能会跟公式对的更整齐一些.

注:

1. 《Reactive Nonholonomic Trajectory Generation via Parametric Optimal Control》此论文中求解的问题Two-
   Point-Boundary-Value problem (TPBVP) , 优化变量为多项式系数, 求解过程很详细值得参考.