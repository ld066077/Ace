---
tags:
  - Readed
---
# 标题::测地声模及残余带状流理论及模拟

## 1.Abstract&Info
### 1.1 Abstract
测地线声模（GAM） 在环向上是对称， 而极向上大致对称。 它是环形等离子体中一种特有的静电扰动模。 测地线声模之所以受到广泛关注， 是因为它可以把漂移波湍流通过非线性相互作用从不稳定的长波区域散射到稳定的短波区域， 进而达到控制和改善等离子体输运的效果。 本文首先调研了测地声模的理论， 包括磁流体力学（MHD） 和动理学理论。 因为托卡马克芯部的温度非常高， 无碰撞的波与粒子相互作用占主导地位， 即非碰撞阻尼起到主要作用， 所以本文采取动理学模拟。 本文使用数值李变换（NLT） 方法， 采取预估校正格式模拟了测地声模随时间的演化。 并且本文采取了改进的场方程求解算法， 可以实现包含磁轴的数值模拟， 以研究测地声模的径向分布问题。 在 RH 测试环节中得到了电场振荡衰减的过程， 也讨论了中心边界条件的选取对结果的影响。 随后研究了在稳定时径向电场的分布， 此过程中修正了初始扰动函数以消除由扰动函数引起的扰动净电荷。 模拟结果显示在扰动区间内与理论仅有比例系数的差别， 在区间外观察到了模拟结果不严格为 0， 这可能是径向传播现象。

### 1.2 Info
**FirstAuthor**:: 肖, 诗麒 
~
**Date**:: 2020
**DOI**: 
**Publication**: 
**PDF**: [肖_测地声模及残余带状流理论及模拟.pdf](file://C:\Users\lyx\Zotero\storage\AUDF5W7E\肖_测地声模及残余带状流理论及模拟.pdf)
**Zotero**: [肖_测地声模及残余带状流理论及模拟.pdf](zotero://select/library/items/AUDF5W7E)


## 2. Annotation%% begin annotations %%


### Imported: 2024-01-13 11:05 晚上


<mark style="background-color: #ff6666">Quote</mark>
>随后研究了在 稳定时径向电场的分布，此过程中修正了初始扰动函数以消除由扰动函数引起 的扰动净电荷。 [jump to](zotero://open-pdf/library/items/AUDF5W7E?page=3&annotation=GFQUVRCT)

标注：修正初始扰动函数怎么做，为什么引起扰动静电荷？

<mark style="background-color: #ff6666">Quote</mark>
>所以带状流可以通过诸如两支波数较高的漂移波产生 一支波数较低的带状流的非线性相互作用获取能量，以抑制漂移波湍流[5-6]。 [jump to](zotero://open-pdf/library/items/AUDF5W7E?page=5&annotation=VQQHQKIH)

标注：非线性相互作用具体是什么

<mark style="background-color: #2ea8e5">Quote</mark>
>一般来说，带状流有两种，低频分支是低频带状流(low frequency zonal flows， 简称 LFZF)[7-8]，或者残余流(residual flows)，稳态带状流(stationary zonal flows)； 有限频率高频分支称为测地线声模(geodesic acoustic mode，简称 GAM)[9-12]。 [jump to](zotero://open-pdf/library/items/AUDF5W7E?page=5&annotation=TPW543PQ)

标注：

<mark style="background-color: #ff6666">Quote</mark>
>托卡马克中某些不稳定性，最典型的是漂移波不稳定性激发的湍流，会激发 m=1，n=0 的密度扰动，产生 m=n=0 的静电扰动。 [jump to](zotero://open-pdf/library/items/AUDF5W7E?page=5&annotation=JRLHCXCF)

标注：漂移波不稳定性激发湍流的图像

<mark style="background-color: #a28ae5">Quote</mark>
>γ绝对值随着 q 的增加先增加后减小，且峰值出现在 q=1.2 处。 因为 q 随着 r 增加而单调增加，所以频率呈内部高外部低，阻尼率中间大两边 小的特点。 [jump to](zotero://open-pdf/library/items/AUDF5W7E?page=10&annotation=LNZHIVLZ)

标注：频率随q递减，阻尼率随q递增后减小

<mark style="background-color: #a28ae5">Quote</mark>
>在内部区域带状流频率为 0，而 q=1.2 之外随着 q 增大频率增大。 阻尼率随着 q 增大而增大。 [jump to](zotero://open-pdf/library/items/AUDF5W7E?page=11&annotation=M3NRVATG)

标注：

<mark style="background-color: #a28ae5">Quote</mark>
>磁流体理论基于单流体方程，对动量方程以磁面平均的方式消除了快磁声波， 并将方程线性化且傅里叶变换，最后得到测地声模和低频带状流的频率解。 [jump to](zotero://open-pdf/library/items/AUDF5W7E?page=11&annotation=8EFWBXP9)

标注：动量方程磁面平均消除快磁声波



%% end annotations %%

## 3.notes
%% begin undefined %%
#### 主要贡献::



%% end undefined %%

%% Import Date: 2024-01-13T23:06:05.866+08:00 %%
