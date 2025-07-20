---
tags:
  - Reading
---
# 标题::A gyrokinetic continuum code based on the numerical Lie transform (NLT) method

## 1.Abstract&Info
### 1.1 Abstract
In this work, we report a novel gyrokinetic simulation method named numerical Lie transform (NLT), which depends on a new physical model derived from the I-transform theory. In this model, the perturbed motion of a particle is decoupled from the unperturbed motion. Due to this property, the unperturbed orbit can be computed in advance and saved as numerical tables for real-time computation. A 4D tensor B-spline interpolation module is developed and applied with the semi-Lagrangian scheme to avoid operator splitting. The NLT code is verified by the Rosenbluth–Hinton test and the linear ITG Cyclone test.

### 1.2 Info
**FirstAuthor**:: Ye, Lei 
**Author**:: Xu, Yingfeng 
**Author**:: Xiao, Xiaotao 
**Author**:: Dai, Zongliang 
**Author**:: Wang, Shaojie 
~
**Date**:: 2016
**DOI**: 10.1016/j.jcp.2016.03.068
**Publication**: Journal of Computational Physics
**PDF**: [Ye et al_2016_A gyrokinetic continuum code based on the numerical Lie transform (NLT) method.pdf](file://E:\Zotero\storage\9V5N5A5T\Ye%20et%20al_2016_A%20gyrokinetic%20continuum%20code%20based%20on%20the%20numerical%20Lie%20transform%20(NLT)%20method.pdf)
**Zotero**: [Ye et al_2016_A gyrokinetic continuum code based on the numerical Lie transform (NLT) method.pdf](zotero://select/library/items/9V5N5A5T)


## 2. Annotation%% begin annotations %%


%% end annotations %%

## 3.notes
%% begin undefined %%
回旋动理学模拟程序 numerical Lie transform（NLT）
通过引入I-变换 求解回旋动理学方程，I-变换本质上是一种特殊的 李变换。
粒子的扰动运动与未扰动运动解耦，由此特性，未扰动运动可以提前的实时计算得到。

漂移波湍流源于温度、密度的空间不均匀性提供的自由能，其决定粒子和热能的传输特性。漂移波湍流是复杂的且本质为非线性，所以解析上非常困难。模拟求解之在实验之外十分有效。

回旋动理学基本思想是解耦粒子快的回旋运动和慢的回旋中心漂移运动[[brizard2007foundations|nonlinearGK]]。对模拟的好处：简化计算，粒子相空间六维->四维回旋中心相空间（though回旋中心分布函数为五维）；由于回旋运动被平均，可增大计算时间步长。


数值求解GKE方法：
1、拉格朗日方法（also PIC）
2、欧拉方法
3、**半拉格朗日方法**（NLT等）：固定的相空间网格点上离散化的分布函数通过每个时间点向后一个时间步长积分扰动轨道，之后插值获得初始点的值来更新。（通常包含扰动轨道积分和相空间插值两部分）

semi-Lagrangian对于回旋中心动力学，跟随完整的回旋中心轨道需要四维插值，使用将高维方程转换为低维方程的 *算子分裂方法* 应用于半拉格朗日法，可以避免高维插值。虽然使得数值更容易实现，但可能引入额外的数值误差，带来杂散耗散（spurious dissipation）。

I-变换的关键思想是将回旋中心的扰动运动与未扰动运动解耦，因此**模拟中只需要未扰动轨道，因为扰动运动已通过 I 变换解耦，我们使用张量积 B 样条进行 4D 插值以避免算子分裂。**

无碰撞等体中，粒子分布函数F(**Z**,t)在回旋中心坐标**Z**=（**X**,$v_\parallel$,$\mu$）满足回旋动理学Vlasov方程
$\frac{\mathrm{d}F}{\mathrm{d}t}=\frac{\partial F}{\partial t}+\dot{\mathbf{X}}\frac{\partial F}{\partial \mathbf{X}}+\dot{v_\parallel}\frac{\partial F}{\partial v_\parallel}=0$
$v_\parallel$为粒子的平行速度，$\mathbf{X}$回旋中心的位置
对时间的全导沿着相空间的回旋中心轨道，因此沿着回旋中心的轨迹分布函数为常量。

未扰动泊松括号，确定未扰动的**导向中心**（GC）的运动。导向中心相空间($\mathbf{X},v_\parallel,\mu,\xi$),$\xi$为回旋角。
这里使用的 I 变换是在短时间间隔 t 内执行的特殊李变换。


%% end undefined %%

%% Import Date: 2024-04-15T16:01:37.921+08:00 %%
