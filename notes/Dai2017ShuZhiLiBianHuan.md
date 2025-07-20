---
tags:
  - Reading
---
# 标题::数值李变换方法及Vlasov系统的非线性数值模拟研究

## 1.Abstract&Info
### 1.1 Abstract
湍流输运是磁约束等离子体研究的重要课题之一。准线性理论很好地预言了一些情况下系统的输运。然而在更多的非线性情况下,由湍流导致的系统输运,在理论上仍然很难处理。随着计算机技术的发展,数值模拟已经在等离子体湍流输运研究中扮演着越来越重要的角色。本文首先介绍了基于全新的数值李变换方法所开发一维无碰撞Vlasov系统的模拟程序。新的程序以I-变换理论为基础,数值上结合了连续性方法以及特征线方法,即避免了粒子模拟存在较大系统噪声的问题,又兼具了特征线方法较强的数值稳定性的特点,在模拟中能采用比较长的时间步长。计算中,我们利用多步变换的方法,在数值上完美地解决了微扰方法在粒子俘获问题上固有的困难,使得程序能够准确地模拟系统线性阶段以及非线性阶段的演化。在一维Landau阻尼问题以及双峰不稳定性问题的算例中,新程序的计算结果与传统方法得到的结果一致。然后我们利用新程序分析了随机电场扰动问题以及双峰不稳定性问题中粒子在速度空间的输运。相较于传统的模拟方法,新的基于微扰理论的模拟方法其中间变量与输运系数有着直接的关联,在数值计算中能够方便地得到系统的输运系数。模拟的结果显示,在系统受到随机场扰动和线性阶段的湍流扰动时,利用新的模拟方法得到的输运系数与实际的结果吻合。而在湍流的非线性阶段,由于大尺度结构的存在,新方法以及准线性方法所计算得到的输运系数与实际的结果有了较大的偏差。接着我们利用数值李变换程序计算了一维无碰撞系统演化过程中的熵产生。与传统数值上用于计算熵的公式不同,我们采用了理论上广为接受的利用粗网平均分布函数来计算系统熵产生的方法。在随机扰动场算例,线性Landau阻尼算例和双峰不稳定性算例中,我们发现随着粗网平均长度的增加,计算所得的熵产生是收敛的,并且当分布函数与麦克斯韦分布接近时,我们计算所得到的熵产生,与热力学定义的熵产生是一致的。我们还讨论了上述情况中粗网平均长度的选取对计算所得熵产生的影响以及非麦克斯韦分布对熵计算的影响。最后,讨论了基于数值李变换方法发展的环位形非线性回旋动理学模拟程序NLT非线性调试中的守恒性问题和滤波问题的处理方法。

### 1.2 Info
**FirstAuthor**:: 戴宗良 
~
**Date**:: 2017
**DOI**: 
**Publication**: 
**PDF**: [戴_2017_数值李变换方法及Vlasov系统的非线性数值模拟研究.pdf](file://C:\Users\lyx\Zotero\storage\UT92D6QT\戴_2017_数值李变换方法及Vlasov系统的非线性数值模拟研究.pdf)
**Zotero**: [戴_2017_数值李变换方法及Vlasov系统的非线性数值模拟研究.pdf](zotero://select/library/items/UT92D6QT)


## 2. Annotation%% begin annotations %%


### Imported: 2024-01-13 10:52 晚上


<mark style="background-color: #2ea8e5">Quote</mark>
>直到目前，准线性理论仍然是理论上计算湍流所导致反常输运的重 要方法。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=17&annotation=F937N59E)

标注：

<mark style="background-color: #2ea8e5">Quote</mark>
>在绝大多数情况下，托卡马克等 离子体湍流是低频的，其特征频率远远低于带电粒子在磁场中的回旋频率。而 系统的宏观输运的时间尺度又远远大于湍流自身的特征时间。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=17&annotation=ZZUMD5DP)

标注：湍流的时间尺度：10^-5s

<mark style="background-color: #2ea8e5">Quote</mark>
>因此在低频湍流的动理学模拟中，计算 粒子的回旋运动是低效且不必要的。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=17&annotation=8FX9NWDR)

标注：

<mark style="background-color: #ff6666">Quote</mark>
>系统的自由度由 整个 6 维相空间缩减到了 4:5 维，同时系统的最大特征频率也大大减小。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=17&annotation=Z2PGYYI9)

标注：4.5维

<mark style="background-color: #2ea8e5">Quote</mark>
>另一类是将系统的 分布函数 f 看作是相空间的连续流体，通过求解 Vlasov 方程模拟系统演化的连 续性方法，其计算结果的噪音远远低于粒子模拟方法所得到的结果。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=18&annotation=Z396IK8K)

标注：

<mark style="background-color: #2ea8e5">Quote</mark>
>连续性方法现有两个分支，一种是直接求解离散的 Vlasov 方程的欧拉方法。 这种方法受限于数值稳定性条件的约束 [38]，计算中的时间步长必须小于空间 网格与速度的比值 ∆t < ∆x/v，而粒子模拟方法并不受限于这个条件。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=18&annotation=DHN6E5F5)

标注：

<mark style="background-color: #2ea8e5">Quote</mark>
>一般情况下，高维插值消耗的时间 非常长，更本无法应用于实际计算。为了克服这一问题，研究者们发展了分裂 方法 [45, 47]，即将高维相空间光滑的特征线，巧妙地分解为由若干个低维相空 间 (一维或二维相空间) 曲线所连接成的折线，这条近似的折线保留了原有曲线 的大部分性质。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=18&annotation=Z5VUFKB3)

标注：

<mark style="background-color: #ff6666">Quote</mark>
>新方法的另一个特点便是依靠未扰动轨道 来计算系统演化，扰动轨道的计算被转化为沿未扰动轨道的积分。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=18&annotation=PZE2AI87)

标注：未扰动轨道积分法与此有什么关系？带电粒子运动的未扰动轨道？

<mark style="background-color: #a28ae5">Quote</mark>
>但是新方法是基于扰动理论得到的，理论上扰动方 法无法处理长时间的非线性共振问题。因此，对于这新方法应用在非线性系统 数值模拟上的具体方法有待研究。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=18&annotation=7RZXLLD2)

标注：

<mark style="background-color: #a28ae5">Quote</mark>
>本文首先发展了数值李变换方法，将 I-变换方法应用于非线性 Vlasov 系统 的模拟。并以此为工具研究了一维 Vlasov 系统的输运性质和无碰撞过程中的熵 增。最后，介绍基于数值李变换方法所发展的托卡马克等离子体非线性回旋动 理学模拟程序 NLT 在非线性调试中的部分问题。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=19&annotation=4TURLGWW)

标注：

<mark style="background-color: #ff6666">Quote</mark>
>(1.4) [jump to](zotero://open-pdf/library/items/UT92D6QT?page=20&annotation=P43JX2NC)

标注：泰勒展开，$\overline{F}(\overline{Z})=F(Z)=T_\epsilon^{-1}F(\overline{Z})$

<mark style="background-color: #ff6666">Quote</mark>
>(1.7) [jump to](zotero://open-pdf/library/items/UT92D6QT?page=20&annotation=MFU5M2LL)

标注：推导过程存疑

<mark style="background-color: #a28ae5">Quote</mark>
>在实际计算中，我们是提 前给定了变换后系统作用量的微分以及哈密顿量的形式，利用方程 (1.11)，从而 确定规范函数以及生成矢量场，进一步确定变换的具体形式。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=21&annotation=8R8H3FJB)

标注：

<mark style="background-color: #2ea8e5">Quote</mark>
>现代回旋动理学通过李变换扰动理论将带电粒子的快变的回旋运动与缓变 的回旋中心的运动解耦。利用同样的思想，I-变换方法在此之上再利用李变换扰 动理论，将回旋中心的运动分解为未扰动运动与扰动引起的偏移 [13]。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=22&annotation=IPQFA68H)

标注：

<mark style="background-color: #ff6666">Quote</mark>
>数值李变换方法使用 I-变换方法的基本方程，而 物理图像却与之不同。正是这点不同使基于微扰方程的数值李变换方法能够模 拟系统的非线性演化。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=25&annotation=PRDLIQE5)

标注：物理图像区别在哪

<mark style="background-color: #a28ae5">Quote</mark>
>假设扰动电场是高频的，由于 离子的惯性远远大于电子，因此这里我们假设其作为背景对扰动没有响应。r [jump to](zotero://open-pdf/library/items/UT92D6QT?page=25&annotation=Y7TD2GBC)

标注：

<mark style="background-color: #ff6666">Quote</mark>
>在多步变换方案中，单次变换所持续的时间 为模拟的时间间隔 τ = ∆t，同时考虑到预估矫正在时间推进上也是二阶精度到， 因此在采用多步变换方案演化系统时，我们不必计算二阶的生成矢量场与规范 函数，但是需要保留一阶生成矢量场所贡献的二阶项。 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=29&annotation=N8VFN4I6)

标注：不必计算二阶生成矢量场和规范函数，但要考虑一阶量贡献的二阶项是否矛盾

<mark style="background-color: #a28ae5">Quote</mark>
>利用数值李变换方法求解 VP 系统分布函数演化的基本流程 [jump to](zotero://open-pdf/library/items/UT92D6QT?page=30&annotation=6HE4BDJY)

标注：



%% end annotations %%

## 3.notes
%% begin undefined %%
#### 主要贡献::



%% end undefined %%

%% Import Date: 2024-01-13T22:52:53.203+08:00 %%
