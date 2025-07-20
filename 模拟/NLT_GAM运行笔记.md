
**maxstep=10** (mesh.dat)   ./input/mesh.dat

**Queue=TH_HPC3N** (sub.sh)-->Queue=debug3N (test)
**Nnodes=16**                          -->      Nnodes=2
### **命令**
[slurm作业调度系统与天河二号的基本命令](https://blog.csdn.net/hong1kun/article/details/120561582)
```
cd XX[/xx] 进入文件
yhq 查看运行序列
./xx.sh  运行sh文件
make NLT 生成.o文件
make clean 删除.o文件
./sub.sh  提交算例
cd ..     返回上一级目录
cd ../..  返回上两级目录
cd或cd ~  返回home目录
run.out   查看运行记录
yhcancel jobid 取消任务
chmod a+x sub.sh   授予权限
cp -r ${HOME}/NLT ./
	The command you provided, `cp -r ${HOME}/NLT ./`, is a Unix/Linux shell command used to copy the directory named "NLT" and its contents from your home directory (`${HOME}`) to the current directory (`./`), recursively (`-r`).

```



---
### notes
$$q(r)=(c_0+c_2*(r/a)^2+...)/\sqrt{1-(r/R_0)^2}$$
```
**读取粒子种类，定义结构类型Par变量ParList，      Global.f90
**PType不等于1 or 2，denn0=1,temt0=1,upaU0=0   Global.f90


**!$OMP PARALLEL DO PRIVATE(i,r)
	- `!$OMP PARALLEL DO`：这是 OpenMP 并行循环的开始标志，它告诉编译器接下来的 `DO` 循环将会并行执行。
	- `PRIVATE(i, r)`：这部分声明了哪些变量是私有的。在这里，`i` 和 `r` 被声明为私有变量，这意味着每个并行线程都有自己的 `i` 和 `r` 变量的副本。
**nofr 径向采样点                              eq_adhoc.in
noftheta
**dr=a/(nofr-1)                              eq_adhoc.f90
r(i) = (i-1)*dr
**F_cn-->q_bar
	if(F_cn == .false.)  then
     do i = 1, nofr
        q_bar(i) = c0 + c2*r(i)**2/a**2
        dpsipdr(i) = Bphi0*r(i)/q_bar(i) !dpsip/dr
        psip(i) = 0.5_dp*a**2*Bphi0/c2*log(1 + c2*r(i)**2/a**2/c0) !psip
     end do
  else
     do i = 1, nofr
        q_bar(i) = c0 + c1*(r(i)/a) + c2*(r(i)/a)**2 + c3*(r(i)/a)**3 + c4*(r(i)/a)**4
        dpsipdr(i) = Bphi0*r(i)/q_bar(i) !dpsip/dr
     end do
     call inrcsnak(r,dpsipdr,cscoef)
     cscoefpsir = cscoef
     psip(1) = 0.0_dp
     do i = 2, nofr
        psip(i) = inrcsitg(r,cscoef,r(1),r(i))
     end do
  end if
**从文件中读取名单中的变量值
	  namelist /global/ R0,a,Bphi0,nofr,noftheta,e0,q0,s0,F_cn,c0,c1,c2,c3,c4, nR, nZ
  !read in the global namelist
  open(10,file='eq_adhoc.in')
  read(10,nml=global)
  close(10)
**问题 iPar=？，parlist(:)   ionNum
**inrcsnak函数给cscoef_q（/g/iota）赋值       eq_adhoc.f90
eqdata.dat 存储q的插值
-->withpsi(3,)
ParNum0=1   <--FlagKE
**`read` 语句会从文件中读取数据，并在读取完成后将文件指针移到下一个位置，准备读取下一行数据
IonNum 离子种类数

**修改particle.dat中PType改变DenN0（仅离子计算需要，电子密度由离子密度计算归一化）计算离子密度方式
**温度和平行速度各自由electron和ion的PType计算得到,TemT0,UpaU0   Global.f90         
**ParList(0)%Ut(i) = ParList(0)%UPar*UpaU0(r,PType) 电子UPar=0
**输出文件：balance_xx.dat为主，（GAM不稳定性关注电场变化，阻尼率和频率）
**`CONTAIN`关键字被用来定义模块的主体部分;`CONTAINS`关键字通常用于定义子程序（subroutine）或函数（function）的实现部分;;`CONTAINS`关键字只能出现在主程序或其他子程序的内部，用于定义嵌套的子程序。在Fortran中，子程序可以被其他子程序调用，有助于模块化和组织代码.
**Equilibrium.f90
	InitEquilibrium
		call LoadEquilA、InitGridA、LoadEquilB
		call InitIonProfile、InitEleProfile、InitPhi0
		call InitPhaseSpace
	InitPhaseSpace
		InitGridB
		RZ_FUNCT
	InitPhi0
	InitEleProfile
	InitIonProfile
	LoadEquilB
	InitGridA
	LoadEquilA
	psi2rx*                   *表示有参数
	rx2psi*
	dpsidrx*
	r2rx*
	rx2r*
	RZ2RxT*
	RxT2RZ*
	FUNCT*
		FT_FUNCT
	RZ_FUNCT
	FT_FUNCT
	get_lap_coef*
		FT_FUNCT
	Phi2coef2d*
	initOrbit
	macroOrbit
	RK4*
	RZRK4*
	dPHIdRZ*
	RZdydt*
	fa_bc*
	fa2flx*
	flx2fa*
**`SUBROUTINE`子例程：子程序没有返回值。它们用于执行一些任务，但不返回结果给调用者。子程序通常通过修改传递给它们的参数来影响调用程序的状态,call语句调用；`FUNCTION`：函数有一个返回值。它们用于执行任务并返回一个值给调用者。函数通常不会修改传递给它们的参数，而是通过 `RETURN` 语句返回一个值
**INTENT
	Fortran 中的 `INTENT` 是用于声明*子程序*（子例程、函数）参数的属性的关键字。它指定了参数在子程序内部的用途，有助于编译器进行编译时的类型检查和优化。
	`INTENT(IN)`：表示参数是输入参数，即在子程序内部只能读取该参数的值，不能修改它。这有助于提高代码的可读性和安全性，因为它明确了参数的作用。
	`INTENT(OUT)`：表示参数是输出参数，即子程序将修改该参数的值并返回给调用者。在调用子程序之前，调用者通常不应该对输出参数赋初值，因为子程序将负责为其赋值。
	`INTENT(INOUT)`：表示参数既是输入参数又是输出参数，子程序可以读取并修改该参数的值。这对于需要在子程序内部修改输入参数的情况很有用。
__dp
	通过添加 `_dp` 后缀，你可以显式地将实数常量声明为双精度，它们通常具有约15位有效数字，提供更高的数值精度。这在科学计算和数值计算中非常有用，特别是当需要更高的数值精度来避免舍入误差和数值不稳定性时
**PhiIN0.dat
*dumpPerTime   mesh.dat 输出时间，normlized by R0/Vth
温度、q  =>  频率、γ


***V4.9
*问题：diagnose 

*field.f90  疑似fieldA and fieldB未使用？
*phi ntheta=16  mesh.dat 极向的采样点？  PhiN0.dat 16相同
	PhiZ.dat 32240/80 -> 403 timelist
	PhiN0.dat/macroPhi.dat  515840/(16*80)->403  timelist
*Timelist.dat  simTime*simTimeFactor
	simTime = simTime + dt  NLT.f90
	InitDt ->dt  Initial.f90
	   dt = dtMax * DToverDTMAX  
	DToverDTMAX = 1.0_dp     and VlasovDF.f90 read DToverDTMAX ?
	dtMax = 2.0 normlized by 2pi/omega_i    4pi
	
	simTimeFactor = 1._dp/R0  Global.f90  ?? 
	L_ref=0.0012
*B_ref=1.5  Bphi0   Global.f90*** and other parameter

新程序V4.9ntheta=16，极向电势并不相同
PhiN0.dat 每同时刻theta向16点径向80点的电势为一组
极向电势*磁面平均*后应该相同，为何不同？

Fortran是**不区分大小写的**
TimeList.dat 不是均匀增加的时间，模拟图像可以正确画出，拟合曲线上纵轴怎么对上？  TIME(i)  403 to 1000

q \times k_r \times rho_i <<1
```


---
## test1-200
### test1 9035536
```
!eq_adhoc.in
&global
        R0 = 1.25
         a = 0.45
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 1.5                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  0.85
        c1 =  0.00
        c2 =  2.41
        c3 =  0.00
        c4 =  0.00
        nR = 513
        nZ = 513
/
```
a.out
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083
```
#### ionnum=2  9064534


---
### test2  9036523
```
!eq_adhoc.in
&global
        R0 = 1.25
         a = 0.45
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 1.5                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1
        c1 =  0.00
        c2 =  3
        c3 =  0.00
        c4 =  0.00
        nR = 513
        nZ = 513
/

```
a.out
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.77905804860493     
  0.890627768262974
```

---

### e_charge=1-test2   9055122
Balance_N00.dat（电子密度剖面）为-1，修改电子charge=1，Balance_N00为1
### ion1 charge=1  test1  ionnum=1   9068441*
**正确的初始分布**
balance_N00.dat因为ion1 charge=-1，所以全为-1
更正后全为1

### test3    9070408
**ionNum=1 去除DenN0 = 1.0_dp，TemT0 = 1.0_dp**
```
!eq_adhoc.in
&global
        R0 = 1.25
         a = 0.45
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 1.5                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  0.85
        c1 =  0.00
        c2 =  2.41
        c3 =  0.00
        c4 =  0.00
        nR = 513
        nZ = 513
/
```
a.out
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083
```
### test4    9071547
**ionNum=1 去除DenN0 = 1.0_dp，TemT0 = 1.0_dp 和 UpaU0=0.0_dp**
```
!eq_adhoc.in
&global
        R0 = 1.25
         a = 0.45
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 1.5                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  0.85
        c1 =  0.00
        c2 =  2.41
        c3 =  0.00
        c4 =  0.00
        nR = 513
        nZ = 513
/
```
a.out
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083
```

### test5    9086113
**ionNum=1 去除DenN0 = 1.0_dp，TemT0 = 1.0_dp 和 UpaU0=0.0_dp**
+
**a=0.60,  R0=1.40**
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.428571428571429     
   1.40000000000000       0.800000000000000        1.48704247446255     
  0.877732473054486
```
### test6      9086233
**ionNum=1 去除DenN0 = 1.0_dp，TemT0 = 1.0_dp 和 UpaU0=0.0_dp**
+
ion charge=2
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083
```

### test7   9086337
**ionNum=1 去除DenN0 = 1.0_dp，TemT0 = 1.0_dp 和 UpaU0=0.0_dp**
+
a=0.75, R0=2 ;c2=3,c3=1,c0=1
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.375000000000000     
   1.40000000000000       0.800000000000000        1.90885428892737     
   1.03643724696304
```

### test8     9086361
**ionNum=1 去除DenN0 = 1.0_dp，TemT0 = 1.0_dp 和 UpaU0=0.0_dp**
+
Bphi0=2
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083
```

### test9   9088613
ionNum=2
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083
```


### test10 
global denn0函数   if(r<0.6a,>0.45a)DenN0=1+delta n

### test11           nnode=16   th   step
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083 

 
### test12       9107363
DenN0=1+delta n
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083
  **maxstep=10**  Nnodes=2  Queue=debug3N (test)

### test13 9107732
16 TH   maxstep=100000
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083 

### test14      9112237
确定r的范围是0.45-0.6
r in0.5-0.55   denn0=1+delta n
nnode=2  debug3N     maxstep=10
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083   
### test15 
r in0.5-0.55   denn0=1+delta n
16 TH   maxstep=100000
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083
```


###  test16 i-edition   
ion  and eletron  PType=1
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083 
```

### test17 i-edition

### test18    9130644
ntheta = 6448  mesh.dat
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083 
```
Segmentation fault

### test19       9131589
r in0.5-0.52   denn0=1+delta n

out：2560 phi

### test20 wed       9132353
r in0.5-0.52   denn0=1+delta n  
q:  0.85  2.41
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083    
```
### test 21 wed   9143270
r in0.5-0.52   denn0=1+delta n   kr=2\*pi/0.02a (698.132)
q: c1-c4=0   c2=0.01  c0 =1.5196881463     eq_adhoc.in
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.54746347038059     
  3.676965618427842E-002

```

### test22 v 4.9  9172476
c0=1.4
无密度扰动
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.42324643888397     
  3.348491112033791E-002

```
### test23  v4.9    9181505
c0=1
Global.f90   DenN0  DenN0 = 1.0_dp+10e-4 sin(r2pi/(0.15))
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.01660459920281     
  3.348491112031839E-002

```

### test24 V4.9   9190381       V4.9原始参数+密度扰动
nrx=80
c0 =  1.40
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        if(r.gt.0.45 .and. r.lt.0.6)then
        DenN0 = 1.0_dp+10e-4sin(r2pi/(0.15))
      else
        DenN0 = 1.0_dp
      end if
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.42324643888397     
  3.348491112033791E-002

```

### test25 V4.9   9190633
c0 =  1.80
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        if(r.gt.0.45 .and. r.lt.0.6)then
        DenN0 = 1.0_dp+10e-4sin(r2pi/(0.15))
      else
        DenN0 = 1.0_dp
      end if

```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.82988827856507     
  3.348491112032405E-002

```
run.out
  InitGlobal time is   3.556323051452637E-002
 psilimitsmall,psilimitbig after normalization=  0.000000000000000E+000
   274.338736366626     
 initPhi are calculated by Balance Equ.
  InitEquilibrium time is    1.67148900032043   end


### test26 V4.9     compared with test21    9194285
q: c1-c4=0   c2=0.01  c0 =1.5196881463     eq_adhoc.in
        if(r.gt.0.45 .and. r.lt.0.6)then
        DenN0 = 1.0_dp+10e-4sin(r2pi/(0.15))
      else
        DenN0 = 1.0_dp
      end if

```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.54746347038059     
  3.676965618427842E-002

```
wrong

### test 27 V4.9  compared with  test20  9196687
q: c1-c4=0   c2=2.41  c0 =0.85     eq_adhoc.in
        if(r.gt.0.45 .and. r.lt.0.6)then
        DenN0 = 1.0_dp+10e-4sin(r2pi/(0.15))
      else
        DenN0 = 1.0_dp
      end if
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.47661818034202     
  0.863089041929083 
```


### 下一步
扫频  Omega随  B   T(0.15-几个kev)和q1-3(或小于1？)的变化
怎么改分布函数？
扰动密度是否应用？ 看code  [[#^2d65a0]]
理论和模拟曲线对比　　Omega　gamma    找模拟/实验的例子
讲故事，来龙去脉  调研GAM  20年的3-5篇理论、实验、模拟文章
调研ETG、EGAM
AI找文章、翻译文章
nrx  80 **2试试？  离子回旋半径量级   模拟结果不变，和理论对比？  [[#^7de7f6]]

### test28 V4.9  9209440

^7de7f6

nrx=160
c0 =  1.40
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        if(r.gt.0.45 .and. r.lt.0.6)then
        DenN0 = 1.0_dp+10e-4sin(r2pi/(0.15))
      else
        DenN0 = 1.0_dp
      end if
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.42324643888397     
  3.348491112033791E-002

```

### test29 V4.9   9209648

^2d65a0

nrx=80
c0 =  1.40
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        if(r.gt.0.45 .and. r.lt.0.6)then
        DenN0 = 1.0_dp+10e-3sin(r2pi/(0.15))
      else
        DenN0 = 1.0_dp
      end if
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.42324643888397     
  3.348491112033791E-002

```

### 勘误 密度扰动1e-4  V4.9  r_max = 0.61   r_min = 0.44


### test30 V4.9   Ti=Te in  0.15Kev-2.25Kev 
nrx=80
c0 =  1.40
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        if(r.gt.0.45 .and. r.lt.0.6)then
        DenN0 = 1.0_dp+1e-4sin(r2pi/(0.15))
      else
        DenN0 = 1.0_dp
      end if
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.360000000000000     
   1.40000000000000       0.800000000000000        1.42324643888397     
  3.348491112033791E-002

```

| 序号 | T_i/kev | job id|
| :-- | :-:|:-:|
|1|0.15|9255475|
|2|0.30|9255408|
|3|0.45|9255644|
|4|0.60|9255657|
|5|0.75|9255666|
|6|0.90| 9257485                |
|7|1.05|       9257674          |
|8|1.20|          9258975      |
|9|1.35|        9259071        |
|10|1.50|      9259201          |
|11|1.65|          9268093     |
|12|1.80|     9268604           |
|13|1.95|         9269866       |
|14|2.10|        9270708        |
|15|2.25|        9270881       |

### test31  newepsi    9285908
c0 =  2.00
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        R0 = 1.71
         a = 0.342
         T_i=T_e=3kev

```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.200000000000000     
   1.40000000000000       0.800000000000000        2.01007563051848     
  1.010101010100414E-002
```

### test32  newepsi    
c0 =  2.00             排除c2设置问题
        c1 =  0.00
        c2 =  0.01
        c3 =  0.00
        c4 =  0.00
        R0 = 1.71
         a = 0.342
         T_i=T_e=3kev

```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.200000000000000     
   1.40000000000000       0.800000000000000        2.01258822505663     
  1.259788900238537E-002

```


fig 1
a 
0.05404	1.69247
0.10814	1.71409
0.16188	1.72484
0.2157	1.73844
0.26944	1.78636
0.32354	1.8305
0.377	1.81129

c
### 下一步
模拟数据处理拟合使用matlab库函数（前一个周期不用，偏后的周期不用）
使用残余带状流的公式作理论图像      ,公式附在图上？
处理数据用PhiZ.dat极向平均后的，PhiN0.dat仅作了环向平均
vlasovDF.f90？**
调 k和调M **
模仿文献作图novikau2017LinearGyrokineticInvestigation （密度扰动仍不变，是否正确？）


### test33      newepsi
nrx=80
c0 =  1.40
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        if(r.gt.0.45 .and. r.lt.0.6)then
        DenN0 = 1.0_dp+1e-4sin(r2pi/(0.15))
      else
        DenN0 = 1.0_dp
      end if
   R0 = 1.71
         a = 0.342
         T=150eV

```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.200000000000000     
   1.40000000000000       0.800000000000000        1.40705294136287     
  1.010101010097337E-002

```



### test34  newepsi
nrx=80
c0 =  1.40
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        if(r.gt.0.45 .and. r.lt.0.6)then
        DenN0 = 1.0_dp+1e-4sin(r2pi/(0.15))
      else
        DenN0 = 1.0_dp
      end if
   R0 = 1.71
         a = 0.342
         T=150eV
  simTotalTime = 80.4
       maxStep = 200000

```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.200000000000000     
   1.40000000000000       0.800000000000000        1.40705294136287     
  1.010101010097337E-002

```

### simtime=35 tmax=87065
### test35 newepsi     q 1-10 1357
a=0.34
R=3.40
epsilon=0.1
Bphi=3.0834
2pi /kr=0.15a
kr rho_i=0.1
T_i=150ev
T_e=3ev   **tau >>me/mi**



### 诊断位置 电场峰点位置(中点)、not密度扰动峰点位置 r=4.8735296459E-0001  使用dat文件中的q和r的数值标注图像

### test36 newepsi   E:\NLTrun\outputV4.9\newepsi vary q
Zotero: [Gao et al_2006_Multiple eigenmodes of geodesic acoustic mode in collisionless plasmas.pdf](zotero://select/library/items/BKZRHJA4)
a=0.34
R=3.40
epsilon=0.1
Bphi0=1.5
2pi/kr=0.15a
kr rho_i=0.1324
T_i=62.1841 ev
T_e=6.21841 ev
mi=2mp
#### try1  9401409
c0=5
c2=0.01
```
 a/R0,q0,s0(input),q0_2,s0_2(output)=  0.100000000000000     
   1.40000000000000       0.800000000000000        5.00876487408918     
  3.505765913999604E-003

```

| 序号 | q(c0= ,c2=0.01)    | job id  |
|----|------|---------|
| 1  | 1.0  |9417192  |
| 2  | 1.5  | 9417634 |
| 3  | 2.0  |  |
| 4  | 2.5  |         |
| 5  | 3.0  |    9421125     |
| 6  | 3.5  |     9421912    |
| 7  | 4.0  |     9422723    |
| 8  | 4.5  |    9422895     |
| 9  | 5.0  |      9423250   |
| 10 | 5.5  | 9423457        |
| 11 | 6.0  |   9424572      |
| 12 | 6.5  | 9425022    |
| 13 | 7.0  |   9425330      |
| 14 | 7.5  |         |
| 15 | 8.0  | 9426852        |
| 16 | 8.5  |         |
| 17 | 9.0  |  9428377       |
| 18 | 9.5  |         |
| 19 | 10.0 |9429790 |
修改径向diag位置，区别不大，0.48->0.525
 
**修改**   newti
kr rho_i=0.100
T_i=35.4987 ev
T_e=3.54987 ev


测试  newti2-4
c0=1.00
1. c2=0    9460426
2. c2=0+T_e=0   9461080
3. T_e=0

#### 修改q剖面为平面

##### test1-0 9466701
T_e=0
q plane
数据不可用，无规振荡
#####  test1  9468321
kr rho_i=0.100
T_i=35.4987 ev
T_e=3.54987 ev
```
q plane
eq_adhoc.f90
qbt = (c0 + c1*(tmpr/a) + c2*(tmpr/a)**2 + c3*(tmpr/a)**3 + c4*(tmpr/a)**4)*sqrt(1-(tmpr/R0)**2)

q_bar(i) = (c0 + c1*(r(i)/a) + c2*(r(i)/a)**2 + c3*(r(i)/a)**3 + c4*(r(i)/a)**4)*sqrt(1-(r(i)/R0)**2)
c0 =  5.00
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00

```
**output**:实频为理论的近似一半，虚频吻合  0.48a处
实频为理论的近似一半，虚频为理论的二倍  0.525a处
##### test2   9468451
T_e=1ev
**output**:实频为理论的近似一半，虚频理论为模拟的1.8倍 0.48a处
实频为理论的近似一半，虚频吻合  0.525a处
##### test3  9468664
```
R_R = 0.55
L_R = 0.45
r_max = 0.56
r_min = 0.44
```
**output**:  实频理论近似为模拟的两倍，虚频理论近似为模拟的3.5倍 0.475a处 
实频为理论的近似一半，虚频吻合  0.5a处
**波的径向传播？实频为何无影响 ->global效应**

##### test4  9473385
lambda=0.1
T_i= 0.01579 KeV
T_e=1eV
**output**:  实频为理论近似1/3,虚频误差25.99%  0.475a处
实频为理论的近似1/3，虚频为理论的近似3/4  0.5a处

##### test5 9473585
lambda=0.1
T_i= 15.79 eV
T_e=0.1 eV
**output**:    0.475a处  拟合出错
实频为理论的近似1/3，虚频为理论的近似32.25%  0.5a处

##### test6 9474153
a=0.68 m
T_i=63.1 eV
T_e=0 eV
B= 1.5 T
k_r rho_i=0.1
**output：** 输出不可用

##### test7  9477020
T_e=1 eV
a=0.68 m
T_i=63.1 eV
B= 1.5 T
k_r rho_i=0.1
**output：** 实频理论为模拟的1.5763倍，虚频理论为模拟的4.4014倍0.475a处
实频与0.475a处几乎不变，虚频为正？且绝对值与理论符合良好 0.5a处
模拟中其他位置处均为实频几乎不变，虚频由负变为正再为负(0.45-0.55)

##### test8
###### try1  9477891   not run
T_i=120 eV
T_e=0 eV
R0 = 2.35
         a = 0.47
         c0 =  2.25
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        Bphi0 = 2.0 
        !Ion1
  mess_i   = 1.0
  k_r rho_i=0.1036
###### try2  9477976  not run
T_i=120 eV
T_e=1 eV
R0 = 2.35
         a = 0.47
         c0 =  2.25
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        Bphi0 = 2.0 
        !Ion1
  mess_i   = 1.0
  k_r rho_i=0.1036
###### try3 9478007 结果竟然还不错
T_i=120 eV
T_e=1 eV
R0 = 2.35
         a = 0.47
         c0 =  5.00
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        Bphi0 = 2.0 
        !Ion1
  mess_i   = 1.0
  k_r rho_i=0.1036
**output：** 实频模拟为理论的1.25倍，虚频模拟为理论的2.22倍 0.5a处
r=0.535a处，振荡先增长又阻尼
r=0.505a处，实频为模拟的1.25倍，虚频模拟为理论的1.5843倍
	进一步缩小画图范围，使得图中仅有一段阻尼率相同的波模，实频模拟为理论的1.2564倍，虚频模拟为理论的1.7931倍
r=0.5a处稍微减小画图区域，使仅包含单一频率，较准确；实频模拟为理论的1.2567倍，虚频模拟为理论的1.6086倍
0.475 a:omg/omgth=1.28715，gam/gamth=2.91946

###### try4 9478112  T_e=0 eV
T_i=120 eV
T_e=0 eV
R0 = 2.35
         a = 0.47
         c0 =  5.00
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        Bphi0 = 2.0 
        !Ion1
  mess_i   = 1.0
  k_r rho_i=0.1036
**output：** 无规则振荡
###### try5  9478512
T_i=120 eV
T_e=1 eV
R0 = 2.35
         a = 0.47
         c0 =  5.00
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        Bphi0 = 2.0 
        !Ion1
  mess_i   = 1.0
  k_r rho_i=0.1036
  r_min=L_R
  r_max=R_R
**output:** 实频模拟为理论的1.25倍，虚频模拟为理论的2.8695倍 0.5a处
###### try6 9481521
T_i=150 eV
T_e=1eV
R0 = 2.35
         a = 0.47
         c0 =  5.00
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        Bphi0 = 2.0 
        !Ion1
  mess_i   = 1.0
  k_r rho_i=
r_min=L_R-0.1
  r_max=R_R+0.1
**output：** 0.475 a:omg/omgth=1.37751，gam/gamth=1.00533
0.5 a:omg/omgth=1.35843，gam/gamth=2.01792
0.525 a:omg/omgth=1.37705，gam/gamth=0.73878
###### scan q
T_i=150 eV
T_e=1eV
R0 = 2.35
lambda_rho=0.1
         a = 0.47
         c0 =  5.00
        c1 =  0.00
        c2 =  0.00
        c3 =  0.00
        c4 =  0.00
        Bphi0 = 2.0 
        !Ion1
  mess_i   = 1.0
  k_r rho_i=
r_min=L_R-0.1
  r_max=R_R+0.1
  
| newT_i | q(c0= ,c2=0.00)    | job id   | ok|
|----|------|---------|-----|
| 1  | 1.0  |9481723  |√|
|2|1.3| 9483521 | √ |
|  |  1.5  | 9481804 |×|
||1.6|  9483044   | × |
| |1.75| 9483443 | × |
||1.84|9483458|×|
| 3  | 2.0  | 9481821 |×|
||2.30| 9483464 |√ |
| 4  | 2.5  |   9481887      |√|
| 5|2.7|9483572 | √|
||2.8|9483547  | × |
|   | 3.0  |     9481913    |×|
| 6  | 3.5  |      9482000   |×|
| | 3.7| 9483632 | √ |
| 7  | 4.0  |     9482024    |√|
| 8  | 4.5  |     9482129    |×|
| |4.7| 9483775  |× |
| |4.8  | 9483790  | √ |
| 9  | 5.0  |         |√|
| 10 |5.3|  9483898  | √ |
| |5.5  |     9482152    |×|
| 11 | 6.0  |    9482163     |×|
|  | 6.3|9485890   | ×  |
| | 6.35   | 9485909  | ×|
| | 6.40|9485927  | √ |
| 12 | 6.5  |   9482220  |√|
| 13 | 7.0  |     9482252    |×|
||7.3 | 9485967  | √ |
| 14 | 7.5  |    9482445     |√|
| 15 | 8.0  |  9482816       |√|
| 16 | 8.5  |  9482874       |√|
| 17 |8.7|9486006  | √ |
| |9.0  |   9482889      |×|
|    | 9.3     |  9485976    |   ×     |
| 18 | 9.5  |  9482903       |√|
| 19 | 10.0 |9483009 |√|

###### dtmax=2.0  -> dtmax=1 运行时间×2
###### try7  9488843
q=4.8
r=0.5a
T_i=150eV
T_e=1eV
Bphi0=2T
dtmax=1.0
**output：** 本来增长的波动，dtmax改为1之后变为阻尼
	0.5 a:omg/omgth=1.37407，gam/gamth=0.678668

###### try8   9490692
if(r.gt.L_R+(R_R-L_R)/4 .and. r.lt.R_R-(R_R-L_R)/4)then
        DenN0 = 1.0_dp+1e-4sin((r-(L_R+(R_R-L_R)/4))2pi2/(R_R-L_R)) 
**output:** 画图拟合区域 : 2-20
	0.5 a:omg/omgth=1.30091，gam/gamth=1.88787
	图像振荡明显仅包含单一模式，时间长后部分长起

###### try9  9491000
if(r.gt.0.49 .and. r.lt.0.51)then
        DenN0 = 1.0_dp+1e-4sin((r-0.49)2pi/(0.02))
**output：** kr rho_i较大 ，0.5
	与理论不符合


###### try10  9491430
扰动波长0.05
if(r.gt.L_R+(R_R-L_R)/4 .and. r.lt.R_R-(R_R-L_R)/4)then
        DenN0 = 1.0_dp+1e-4sin((r-(L_R+(R_R-L_R)/4))2pi2/(R_R-L_R)) 
T_i=29eV
T_e=1eV
**output:** 结果不合理，增长
###### try11  9491780
a=0.235
R0=1.175
T_i=120eV
T_e=1eV
B_phi=8T
扰动波长=0.05a
**output:** 结果不合理，增长

###### try12   9492191
VlasovDF.f90
GAM_B Init
FALSE -> TRUE
Random Init
TRUE -> FALSE
扰动波长 0.1
dtmax=2.0
scan q参数



### 目标：模拟符合理论，在0.5a处，kr rho_i较小；not run 卡在equilibrium时间->调整q   ；


### 拟合不准确应该调整left_dot=2 ref_dot=3组合; r/R(0.1) not eq a/R(0.2)
### T_e为0eV出现无规则振荡

### test37 验证理论符合 yanzheng_newepsi  9422848
Zotero: [Zhao et al_2019_Linear gyrokinetic simulations of zonal flows in toroidal rotating plasmas.pdf](zotero://select/library/items/HQFH44EB)
nrx=72 ntheta=64  nvpara=156   nmu=64
R_R = 0.55   L_R = 0.45  r_min=0.44  r_max=0.56                               mesh.dat
R0 = 2.35   a = 0.47
Bphi0 = 2.0
c0 =  2.25
        c1 =  0.00
        c2 =  0.001
        c3 =  0.00
        c4 =  0.00           eq_adhoc.in                    
mess_i=1 
T_i=120ev 
T_e=0 ev
kr rho_i=0.1058                particle.dat
#### grid origin ed  9427125   9430570S
not run

#### q flat  网格设置
##### not run
##### r_min=L_R  r_max=R_R 9476422
not run

### Matlab删除数组数据
```
% 创建示例数组 
data = [12, 45, NaN, 23, 56, NaN, 34, 67]; % 检测数组中是否存在NaN值 
nan_indices = isnan(data);
% 删除包含NaN的数据 
data = data(~nan_indices);
% 显示处理后的数组
disp(data);
```


### test38   modified wed 

output: 0.5 a:omg/omgth=1.24606，gam/gamth=3.38061

### test39 标准重复 9493535  wrong r/R
a/R=0.1   a=0.5 R=5.0  × 
Bphi0=2
nrx=60
C0=2.25
q有梯度
Te=0.001  Ti=0.15
rmin=0.69 rmax=0.81
RL=0.70   RR=0.80
DenN0= (r-0.70)/0.1   0.70-0.80
diagr=0.72
re=2.6084 im=-0.1414
simulation:re=2.487 im=-0.1461 
mess_i=2
k_r rho_i=0.1573
**output：** 0.72 a:omg/omgth=0.979533，gam/gamth=0.612122
0.725 a:omg/omgth=0.976775，gam/gamth=0.577991
0.75 a:omg/omgth=0.992149，gam/gamth=-0.0512746
0.775 a:omg/omgth=0.977615，gam/gamth=0.586447
### test40 biaozhun   9493816  wrong r/R
a/R=0.1   a=0.5 R=5.0×
Bphi0=2
nrx=60
C0=2.25
q有梯度
Te=0.001  Ti=0.15
rmin=0.69 rmax=0.81
RL=0.70   RR=0.80
DenN0= (r-0.70)/0.1   0.70-0.80
diagr=0.72
re=2.6084 im=-0.1414
simulation:re=2.487 im=-0.1461 
mess_i=1
k_r rho_i=0.1112
**output：** 

### test41 biaozhun   9493989
mess_i=1
**output：** 0.72 a:omg/omgth=1.4067，gam/gamth=1.36725
	0.75 a:omg/omgth=1.38487，gam/gamth=2.10632
	0.775 a:omg/omgth=1.37657，gam/gamth=0.715881
	0.725 a:omg/omgth=1.39469，gam/gamth=0.0965115

### test42 biaozhun 9494128
mess_i=2
**output：** **0.75 a:omg/omgth=0.967131,gam/gamth=0.841044**
	0.775 a:omg/omgth=0.984566，gam/gamth=0.292533
	0.72 a:omg/omgth=0.984516，gam/gamth=0.270205
	0.725 a:omg/omgth=0.984574，gam/gamth=0.393091

### test43 biaozhun  9494511
mess_i=2
VlasovDF.f90
GAM_B Init
FALSE -> TRUE
Random Init
TRUE -> FALSE
**output：** 0.725 a:omg/omgth=0.980867 gam/gamth=0.286028
	0.775 a:omg/omgth=0.976769，gam/gamth=0.14664
	0.72 a:omg/omgth=0.987437，gam/gamth=-0.0844087
	**0.75 a:omg/omgth=0.968743，gam/gamth=1.03887**

### test44 biaozhun 9494870
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
q 有梯度，c0=2.25
**output** not run

0.75a/R != 0.5a/R ? ->test47
### test45 biaozhun 9496294
mess_i=2
q有梯度，c0=2.25×  3×  c0=4  √
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           R0=5.0
           a=0.5
**output：** 0.5 a:omg/omgth=0.985756，gam/gamth=1.65321


### test46 biaozhun 9497247
q flat，c0=2.25×   3 ×  2.2563 ×    4√
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
**output:** 0.5 a:omg/omgth=0.984416，gam/gamth=1.66899

### test47 biaozhun 9496891

^9b9386

mess_i=2
q有梯度，  c0=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
**output: 完美**0.5 a:omg/omgth=0.968135，gam/gamth=1.04997


### test48 biaozhun  9497506
q 有梯度
c0=3
**output：** 0.5 a:omg/omgth=0.985777，gam/gamth=1.23019

### test49 biaozhun 重复[test47](#^9b9386)   9498077
一模一样

### test50 biaozhun  failed
q   flat, scan q from 1 to 10
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           mess_i=2


### test51 biaozhun 9503701
mess_i=2
q有梯度，c0= 9.5
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           R0 = 3.34
           a=0.5
**output：** 0.5 a:omg/omgth=0.971763，gam/gamth=1.79119

### test52 biaozhun/modi ed(Debugq)  9513423
compared with 45/46
mess_i=2
q有梯度，c0= 4
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           R0 = 3.34
           a=0.5
**output：** 0.5 a:omg/omgth=0.98361，gam/gamth=1.6721

### test53 debug（Debugq） 9512815
print *, ' begin1'

    allocate(h_Lambda_alpha_c(nrx))

    !$OMP PARALLEL DO PRIVATE(i,m,n,tmp)

    do i = 1, nrx

      n = (nFltr-1) * ntm

      tmp = 0._dp

      do while (tmp .lt. mnFactorC)

        n = n - ntm

        do m = NINT(n*mesh_q(i)-dmFltr), NINT(n*mesh_q(i)+dmFltr)

          call comput_filter_factor(m,n,i,tmp)

          if (tmp .ge. mnFactorC) exit

          print *, ' tmp=',tmp

        end do

      end do

      n = n/ntm

      h_Lambda_alpha_c(i) = nkc*real(n,dp)/(lalph*FactorRe)

    end do

    !$OMP END PARALLEL DO

    print *, ' after1'
修改为->
!!$OMP PARALLEL DO PRIVATE(i,m,n,tmp)
    do i = 1, nrx
    !  n = (nFltr-1) * ntm
     ! tmp = 0._dp
     ! do while (tmp .lt. mnFactorC)
      !  n = n - ntm
       ! do m = NINT(n*mesh_q(i)-dmFltr), NINT(n*mesh_q(i)+dmFltr)
        !  call comput_filter_factor(m,n,i,tmp)
         ! if (tmp .ge. mnFactorC) exit
        !end do
      !end do
     ! n = n/ntm
      h_Lambda_alpha_c(i) = 10
    end do
    !!$OMP END PARALLEL DO


### test54 modi ed  
q   flat, scan q from 1 to 5
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           mess_i=2
           T_i=150
           T_e=1
           nrx=60
[output](file:///F:%5CSimulation%5Cbiaozhun%5Cmodi%20ed%5Ctest54%5Coutput)
scanq
kr rho_i = 0.157236

| 序号(文件夹) | q(flat) | job id  |
|---------|---------|---------|
| 1       | 1       | 9513551 |
| 2       | 1.25    | 9513661 |
| 3       | 1.5     | 9513788 |
| 4       | 1.75    | 9513862 |
| 5       | 2       | 9513956 |
| 6       | 2.25    | 9515524 |
| 7       | 2.5     | 9515571 |
| 8       | 2.75    | 9515629 |
| 9       | 3       | 9515681 |
| 10      | 3.25    | 9515711 |
| 11      | 3.5     | 9515768 |
| 12      | 3.75    | 9515848 |
| 13      | 4       | 9516310 |
| 14      | 4.25    | 9516387 |
| 15      | 4.5     | 9516404 |
| 16      | 4.75    | 9516691 |
| 17      | 5       | 9516775 |

### test55  modi ed
q flat
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           mess_i=2
ntheta=32

2.25  9518955
2.50  9519056
2.75  9522697
3.00  9522826

### test56 modi ed

ntheta=64

2.25  9527546
2.50  9527120
2.75  9526854
3.00  9524617

### test57  modi ed  9527936
T_e=0
ntheta=16
c0=1.75
**output** 高频振荡

### test58 modi ed 9530361
q   flat ,q=1.75
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=150
           T_e=1
a* 0.12/rho_i=67.8071
**output** 0.5 a:omg/omgth=1.41267，gam/gamth=1.45488
### test59 modi ed
$k_r\rho_i=0.0524$
q   flat ,q=1--6
r_max = 0.66
         r_min = 0.34
         R_R = 0.65
           L_R = 0.35
           a=0.5
           R0=3.34
           **mess_i=2**
         T_i=150
           T_e=1
        nrx=130

| 序号(文件夹) | q(flat) | job id  |
|---------|---------|---------|
| 1       | 1       | 9532415 |
| 2       | 1.25    | 9552517 |
| 3       | 1.5     | 9552311 |
| 4       | 1.75    | 9530781 |
| 5       | 2       | 9552163 |
| 6       | 2.25    | 9551016 |
| 7       | 2.5     | 9545404 |
| 8       | 2.75    | 9544130 |
| 9       | 3       | 9542449 |
| 10      | 3.25    | 9542276 |
| 11      | 3.5     | 9542114 |
| 12      | 3.75    | 9541858 |
| 13      | 4       | 9541651 |
| 14      | 4.25    | 9541169 |
| 15      | 4.5     | 9539236 |
| 16      | 4.75    | 9534988 |
| 17      | 5       | 9533672 |
| 18      | 5.25    | 9532969 |
| 19      | 5.5     | 9532879 |
| 20      | 5.75    | 9532826 |
| 21      | 6       | 9532570 |


### test60 modi ed  9554273
$k_r\rho_i=0.1572$
q   flat ,q=1.75
r_max = 0.56
         r_min = 0.44
	         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **$mess_i=2$**
         T_i=150
           T_e=1
        nrx=40
        B_phi0=2.0
ionnum=2
&Par 1&2
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.15
**output** simulation：$\gamma=-0.113868C_s/R_0\quad\omega=2.17941C_s/R_0$
theory：$\gamma=-0.0995763C_s/R_0\quad\omega=2.19594C_s/R_0$
### test61 modi ed  9555210

^79a0f3

**Den=0.5**
	$k_r\rho_i=0.1572$
q   flat ,q=1.75
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=2**
         T_i=150
           T_e=1
        nrx=40
        B_phi0=2.0
ionnum=2
&Par 1&2
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 0.5**
  Upa    = 0.0
  Tem    = 0.15
**output** simulation：$\gamma=-0.108768C_s/R_0\quad\omega=2.17726C_s/R_0$
对比单个粒子模拟结果:$\gamma=-0.108458\quad\omega=2.18024$
theory：$\gamma=-0.0995763C_s/R_0\quad\omega=2.19594C_s/R_0$
### test62 modi ed 9559012
$k_r\rho_i=0.1112$
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=150
           T_e=1
        nrx=60
        B_phi0=2.0
ionnum=1
&Par 1
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 1.0**
  Upa    = 0.0
  Tem    = 0.15
**output** 0.5 a:omg/omgth=1.38108，gam/gamth=1.53294

### test63 modi ed 9559709
$k_r\rho_i=0.1112$
q   flat ,q=1.75
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=150
           T_e=1
        nrx=60
        B_phi0=2.0
ionnum=1
&Par 1
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 1.0**
  Upa    = 0.0
  Tem    = 0.15
**output** 
0.5 a:omg/omgth=1.41578，gam/gamth=1.50129
### test64 modi ed 9559820
$k_r\rho_i=0.1112$
q   flat ,q=3.75
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=150
           T_e=1
        nrx=60
        B_phi0=2.0
ionnum=1
&Par 1
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 1.0**
  Upa    = 0.0
  Tem    = 0.15


### test65 modi ed 9560598

$k_r\rho_i=0.1284$
q   flat ,q=3.75
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=200
           T_e=1
        nrx=60
        B_phi0=2.0
ionnum=1
&Par 1
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 1.0**
  Upa    = 0.0
  Tem    = 0.20
  **output** 0.5 a:omg/omgth=1.61434，gam/gamth=2.29244
### test66 modi ed 9561190

^54f934

$k_r\rho_i=0.0995$
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=120
           T_e=1
        nrx=60
        B_phi0=2.0
ionnum=1
&Par 1
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 1.0**
  Upa    = 0.0
  Tem    = 0.12
**output** 0.5 a:omg/omgth=1.23806，gam/gamth=1.36075
### test67 modi ed 9561266
$k_r\rho_i=0.0908$
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=100
           T_e=1
        nrx=60
        B_phi0=2.0
ionnum=1
&Par 1
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 1.0**
  Upa    = 0.0
  Tem    = 0.10
  **output** 0.5 a:omg/omgth=1.13192，gam/gamth=1.22415

### 为什么换了离子，温度也要换才能与理论符合？？

### test68 modi ed 9561348
$k_r\rho_i=0.0861$
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=90
           T_e=1
        nrx=60
        B_phi0=2.0
ionnum=1
&Par 1
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 1.0**
  Upa    = 0.0
  Tem    = 0.09
**output** 0.5 a:omg/omgth=1.07461，gam/gamth=1.14986

### test69 modi ed 9564521

^1f3025

$k_r\rho_i=0.0812$
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
ionnum=1
&Par 1
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 1.0**
  Upa    = 0.0
  Tem    = 0.08
**output完美** 0.5 a:omg/omgth=1.01483，gam/gamth=1.09759

### 改rho后修改nrx

### test70  compared with [test66](#^54f934)  9563710
$k_r\rho_i=0.0995$
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=120
           T_e=1
        nrx=80
        B_phi0=2.0
ionnum=1
&Par 1
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 1.0**
  Upa    = 0.0
  Tem    = 0.12
  **output** 0.5 a:omg/omgth=1.23881，gam/gamth=1.40577

### test71 加杂质from[test69](#^1f3025)  modied H(He 10%)9566338
$k_r\rho_i=0.0812$!!!!!
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
        H: Den=0.818
ionnum=2
&Par 2
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  **Den    = 0.0909**
  Upa    = 0.0
  Tem    = 0.08

### 有效电荷 电荷密度比
$f_{Cs}=n_sZ_s/n_e$
$\sum_s f_{Cs}=1$
$Z_{eff}=\sum_{s\neq e}f_{Cs}Z_s$

### test72 time * 3  9568430 modi ed
simTotalTime = 40.2  ->  120.6
       maxStep = 100000  -> 300000
$k_r\rho_i=0.0812$!!!!!
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
        H: Den=0.818
ionnum=2
&Par 2
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  **Den    = 0.0909**
  Upa    = 0.0
  Tem    = 0.08

### test73 time * 2 modi ed  9571064
simTotalTime = 40.2  ->  80.4
       maxStep = 100000  -> 200000
$k_r\rho_i=0.0812$
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
ionnum=1
&Par 1
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 1.0**
  Upa    = 0.0
  Tem    = 0.08
**output** 0.5 a:omg/omgth=1.01502，gam/gamth=1.11394

### test74 modi ed  9575196

^74a1df

simTotalTime = 40.2  ->  80.4
       maxStep = 100000  -> 200000
$k_r\rho_i=0.0812$!!!!!
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
        H: Den=0.818
ionnum=2
&Par 2
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  **Den    = 0.0909**
  Upa    = 0.0
  Tem    = 0.08

### test75 modi ed  $T_{He}/T_{H}=10$  9575998

^4c20c3

simTotalTime = 40.2  ->  80.4
       maxStep = 100000  -> 200000
$k_r\rho_i=0.0812$!!!!!
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
        H: Den=0.818
ionnum=2
&Par 2
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  **Den    = 0.0909**
  Upa    = 0.0
  Tem    = 0.80
**output** GAM 在21$R/C_s$内衰减

### test76 Ar(0.6%) modi ed 9580620 not run

^05a41e

Ionnum = 3
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.001
/

&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.818
  Upa    = 0.0
  Tem    = 0.08
/

&Par
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.0909
  Upa    = 0.0
  Tem    = 0.80
/
%默认为完全电离
电荷数满足,自动满足 $\sum_sn_s*e_s=-n_e*e$
&Par
  !Ion3
  mess   = 40.0
  charge = 0.0 **错**
  vMax   = 4.0
  PType  = 4
  Den    = 0.006
  Upa    = 0.0
  Tem    = 0.08
/

### test77  modi ed  $T_{He}/T_{H}=5$  9583097

^dbde54

simTotalTime = 40.2  ->  80.4
       maxStep = 100000  -> 200000
$k_r\rho_i=0.0812$!!!!!
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
        H: Den=0.818
ionnum=2
&Par 2
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  **Den    = 0.0909**
  Upa    = 0.0
  Tem    = 0.40

### test78 Ar(0.28%) modi ed 9583715
Ionnum = 3
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.001
/

&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.818
  Upa    = 0.0
  Tem    = 0.08
/

&Par
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.0909
  Upa    = 0.0
  Tem    = 0.80
/

电荷数满足 $\sum_sn_s*e_s=-n_e*e$
&Par
  !Ion3
  mess   = 40.0
  charge = 0.0  **错**
  vMax   = 4.0
  PType  = 4
  Den    = 0.0028
  Upa    = 0.0
  Tem    = 0.08
/

### test79 modi ed compared with [test61](#^79a0f3) 9587021

^0c88e4

VlasovDF.f90
	!!!!!!!!!!!!!!!!!!

          !!  GAM_B Init  !!

          !!!!!!!!!!!!!!!!!!

          if (.TRUE.) then

            ! initial perterbation for RH test

            !if (iPar .eq. 1) then

              !$OMP PARALLEL DO PRIVATE(i,j,k,tmpr)

              do i = 1, nrx

                tmpr = rx2r(Vrx(i))

                if ( (tmpr .ge. L_R) .and. (tmpr .le. R_R) ) then

                  do k = 1, ntheta

                  do j = 1, nalpha

                    Package4MPI(i,j,k) = sin((tmpr-L_R)*twopi/(R_R-L_R))/Jcb_FA(i)

                  end do

                  end do

                end if

              end do

              !$OMP END PARALLEL DO

            !end if

          end if

**Den=0.5**
	$k_r\rho_i=0.1572$
q   flat ,q=1.75
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=2**
         T_i=150
           T_e=1
        nrx=40
        B_phi0=2.0
ionnum=2
&Par 1&2
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 0.5**
  Upa    = 0.0
  Tem    = 0.15
  **output** simulation：$\gamma=-0.113868C_s/R_0\quad\omega=2.17941C_s/R_0$
theory：$\gamma=-0.0995763C_s/R_0\quad\omega=2.19594C_s/R_0$

### 同时提交多个任务，run后再提交修改文件，或者建立多个文件夹

### test80 modi ed compared with [test79](#^0c88e4) 9587171 9587406
**Den=1.0**
	$k_r\rho_i=0.1572$
q   flat ,q=1.75
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=2**
         T_i=150
           T_e=1
        nrx=40
        B_phi0=2.0
ionnum=1
&Par 1
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  **Den    = 1.0**
  Upa    = 0.0
  Tem    = 0.15
  **output** simulation：$\gamma=-0.113868C_s/R_0\quad\omega=2.17941C_s/R_0$
theory：$\gamma=-0.0995763C_s/R_0\quad\omega=2.19594C_s/R_0$
对比说明之前多粒子的密度扰动没有有效得到计算
### test81 modi ed compare with [test74](#^74a1df) 9587554

^29c0dc

simTotalTime = 40.2 
       maxStep = 100000 
$k_r\rho_i=0.0812$!!!!!
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
        H: Den=0.818
ionnum=2
&Par 2
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  **Den    = 0.0909**
  Upa    = 0.0
  Tem    = 0.08

### test82 modi ed [test75](#^4c20c3) 9587585
simTotalTime = 40.2 
       maxStep = 100000 
$k_r\rho_i=0.0812$!!!!!
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
        H: Den=0.818
ionnum=2
&Par 2
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  **Den    = 0.0909**
  Upa    = 0.0
  Tem    = 0.80

### test83 modi ed [test77](#^dbde54) 9587602
simTotalTime = 40.2 
       maxStep = 100000
$k_r\rho_i=0.0812$!!!!!
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
        H: Den=0.818
ionnum=2
&Par 2
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  **Den    = 0.0909**
  Upa    = 0.0
  Tem    = 0.40

### test84 modi ed Ar(0.28%) He(10%) 9589169
$T_{He}/T_{H}=5$  
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.001
/

&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.8972
  Upa    = 0.0
  Tem    = 0.08
/

&Par
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.1
  Upa    = 0.0
  Tem    = 0.40
/

&Par
  !Ion3
  mess   = 40.0
  charge = 8.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.0028
  Upa    = 0.0
  Tem    = 0.08
/

### test85 modi ed Ar(0.6%) He(10%) 9589965
$T_{He}/T_{H}=5$  
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.001
/

&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.8940
  Upa    = 0.0
  Tem    = 0.08
/

&Par
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.1
  Upa    = 0.0
  Tem    = 0.40
/

&Par
  !Ion3
  mess   = 40.0
  charge = 8.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.006
  Upa    = 0.0
  Tem    = 0.08
/

### test86 modi ed Ar(0.28%) He(10%)9591204
$T_{He}/T_{H}=10$
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.001
/

&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.8972
  Upa    = 0.0
  Tem    = 0.08
/

&Par
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.1
  Upa    = 0.0
  Tem    = 0.80
/

&Par
  !Ion3
  mess   = 40.0
  charge = 8.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.0028
  Upa    = 0.0
  Tem    = 0.08
/

### test87 modi ed Ar(0.6%) He(10%) 9591220
$T_{He}/T_{H}=10$
**三种离子密度设置**：Ar：$n_{Ar}/\sum_sn_s=0.6\%$
He：$n_{Ar}/\sum_sn_s=10\%$
$n_e$设为1，但是运行结果不为1，根据电荷数自动满足 **存疑**
**因此$n_e$准确应该改为$\sum_sn_s$**
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.001
/

&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.8940
  Upa    = 0.0
  Tem    = 0.08
/

&Par
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.1
  Upa    = 0.0
  Tem    = 0.80
/

&Par
  !Ion3
  mess   = 40.0
  charge = 8.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.006
  Upa    = 0.0
  Tem    = 0.08
/
**output** Ar浓度增加后，阻尼率增大，频率减小，residual flow稍有抬升

|   |   |   |   |
|---|---|---|---|
||residual flow|gamma|omega|
|H|0.029705|-0.03655|2.09672|
|H(He0.1)|0.029799|-0.02386|1.900543|
|H(He0.1)T_He/T_H=5|0.03335|-0.08482|2.147078|
|H(He0.1)T_He/T_H=10|0.034866|-0.25209|2.486466|
|H(He0.1)T_He/T_H=10 Ar=0.28%|0.035936|-0.25536|2.406883|
|H(He0.1)T_He/T_H=10 Ar=0.6%|0.036985|-0.26259|2.326357|

### test88 modi ed 9593228  [test81](#^29c0dc)
FlagKE = .F. -> .T.
dtMax = 2.0-> 0.1
simTotalTime = 40.2 
       maxStep = 100000 
$k_r\rho_i=0.0812$!!!!!
q   flat ,q=2.25
r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
        H: Den=0.818
ionnum=2
&Par 2
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 4
  **Den    = 0.0909**
  Upa    = 0.0
  Tem    = 0.08

### 计算密度 修正
密度比有影响，密度整体大小无影响，但密度比为与电子的密度比
n_He=0.1 ne
nH=ne-2 n_He

### test89 modi ed 9623135
q   flat ,q=2.25
simTotalTime = 40.2 
       maxStep = 100000 
       r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
        H: Den=0.8
        ionnum=2
&Par 2
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 5
  **Den    = 0.1**
  Upa    = 0.0
  Tem    = 0.08

```外部密度高，内部低
if (PType .eq. 5)then

      rm     = 0.40_dp

      deltaN = 0.05_dp

      kappaN = 2.2320_dp

      DenN0 = exp(kappaN*ra/R0*( r - rm &

                   - deltaN*tanh((r-Vmr(1)  )/deltaN) &

                   - deltaN*tanh((r-Vmr(nrx))/deltaN) ))

    end if
```
### test90 modi ed 9627160
q   flat ,q=2.25
simTotalTime = 40.2 
       maxStep = 100000 
       r_max = 0.56
         r_min = 0.44
         R_R = 0.55
           L_R = 0.45
           a=0.5
           R0=3.34
           **mess_i=1**
         T_i=80
           T_e=1
        nrx=80
        B_phi0=2.0
        H: Den=0.8
        ionnum=2
&Par 2
  !Ion2
  mess   = 4.0
  charge = 2.0
  vMax   = 4.0
  PType  = 5
  **Den    = 0.1**
  Upa    = 0.0
  Tem    = 0.08
```杂质聚芯
  if (PType .eq. 5)then

      rm     = 0.40_dp

      deltaN = 0.05_dp

      kappaN = 2.2320_dp

      DenN0 = exp(-kappaN*ra/R0*( r - rm &

                   - deltaN*tanh((r-Vmr(1)  )/deltaN) &

                   - deltaN*tanh((r-Vmr(nrx))/deltaN) ))

    end if
```

### code order
#### 代码编号
Debugq  9  modi ed
edgeimpurity 10
rotation 11
EPs1.0  12 ^fd13c9
SDFtest 13  slowing+double
SDFtest ->SDF1(impurity queue)   SDF2(double shifted queue)  SDF3(slowing down queue)
#### 代码迭代关系
edgeimpurity 10<- Debugq 9
edgeimpurity 10-> EPs1.0
EPs1.0-->SDFtest

### test91 10 9632734
      FlagAxis = .F.
         r_max = 0.91
         r_min = 0.79
        rxType = 2
           nrx = 80
           ntm = 1
        nalpha = 2
        ntheta = 16
        nvpara = 64
        muType = 3
           nmu = 16
         dtMax = 2.0
  simTotalTime = 40.2
       maxStep = 100000
           R_R = 0.90
           L_R = 0.80
    IonNum = 1
      !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.08

**output** 0.85 a:omg/omgth=1.00831，gam/gamth=1.06433

### test92 10 9632799
  Tem    = 0.10

### test93 10 9635519
a=0.5
R0=4.25
  Tem    = 0.10

### test94 10 9636732
q=3
a=0.5
R0=4.25
  Tem    = 0.10

### test95 11 9638048
UOMEG=0.6018

### test96 10 9638121 

^ade144

q=3
a=0.5
R0=4.25
Tem = 0.08
**output** omega =2.03 
0.85 a:omg/omgth=1.01873，gam/gamth=1.35407
实频符合Ye2013(11a)，虚频不符合

### test97 11 9638224
UOMEGA=0.1

### test98 10 9638238
q=6
a=0.5
R0=4.25
Tem = 0.08
### test99 11 9638261
UOMEGA=0.0

### test100 10 9638383
q=5
a=0.5
R0=4.25
Tem = 0.06
**output**  0.85 a:omg/omgth=0.88411，gam/gamth=1.74547
### test101 10 9638570
q=2.25
a=0.5
R0=4.25
Tem = 0.06

### test102 10 9640784
q=3
a=0.5
R0=4.25
Tem = 0.06

### test103 10 9640916
q=4
a=0.5
R0=4.25
Tem = 0.06

### test104 10 9641111
q=4
a=0.4
R0=3.25
Tem = 0.06

### test105 10 9641398
a=0.5
R0=3.34
Tem = 0.06
q=4
**output** 0.85 a:omg/omgth=0.882447，gam/gamth=0.889502
### test106 10 9641942
a=0.3
R0=3.34
Tem = 0.06
q=4

### test107 10 9645998
q=2.25
a=0.5
R0=4.25
Tem = 0.06

### test108 10 9647440
q=3
a=0.5
R0=4.25
Tem = 0.06

### test109 10 [test96](#^ade144) 9649347
q=3
a=0.5
R0=4.25
Tem = 2.00
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 2.00

### test110 10 9650320
q=3.5
B0=8.2T
R0=2.35m
a=0.47m
T_i=2keV
nrx = 80
           ntm = 1
        nalpha = 2
        ntheta = 16
        nvpara = 64
        muType = 3
           nmu = 16
!nrx=64
!ntheta=128
!nvpara=128
!nmu=32

### test111 10 9650748
q=3.5
B0=8.2T
R0=2.35m
a=0.47m
T_i=2keV
nrx=64
ntheta=128
nvpara=128
nmu=32

### test112 10 9653309
q=3.5
B0=8.2T
R0=2.35m
a=0.47m
T_i=0.8keV
nrx=80
ntheta=16
nvpara=64
nmu=16

### test113 10 9653826
q=3.5
B0=2.0T
R0=2.35m
a=0.47m
T_i=0.8keV
nrx=80
ntheta=16
nvpara=64
nmu=16

### test114 10 9653996
q=3.5
B0=8.2T
R0=2.35m
a=0.47m
T_i=0.8keV
nrx=100
ntheta=16
nvpara=64
nmu=16

### test115 10 9654362
q=3.5
B0=8.2T
R0=2.35m
a=0.47m
T_i=0.2keV
nrx=200
ntheta=16
nvpara=64
nmu=16

### test116 10 9655795
q=3.5
B0=2T
R0=2.35m
a=0.47m
T_i=0.2keV
nrx=80
ntheta=16
nvpara=64
nmu=16

### test117 10 9655989
q=3.5
B0=2T
R0=2.35m
a=0.47m
T_i=0.1keV
nrx=80
ntheta=16
nvpara=64
nmu=16

### test118 10  9656551
q=3.5
B0=2T
R0=2.35m
a=0.47m
T_i=0.15keV
nrx=80
ntheta=16
nvpara=64
nmu=16

### test119 10 9657567
q=2.25
B0=2T
R0=2.35m
a=0.47m
T_i=0.15keV
nrx=80
ntheta=16
nvpara=64
nmu=16

### test120 10 9660245
q=1.5
B0=2T
R0=2.35m
a=0.47m
T_i=0.15keV
nrx=80
ntheta=16
nvpara=64
nmu=16

### test121 10 9660842
q=5
B0=2T
R0=2.35m
a=0.47m
T_i=0.15keV
nrx=80
ntheta=16
nvpara=64
nmu=16

### N0处不需加delta n扰动->VlasovDF

### test122 10 9662935
q=5
B0=2T
R0=2.35m
a=0.47m
T_i=0.15keV
nrx=80
ntheta=16
nvpara=64
nmu=16

### test123 10 9665515
q=3
B0=2T
R0=2.35m
a=0.47m
T_i=0.15keV
nrx=80
ntheta=16
nvpara=64
nmu=16

### test124 10 9666949
q=3
B0=2T
R0=2.35m
a=0.47m
T_i=0.1keV
nrx=80
ntheta=16
nvpara=64
nmu=16
**output** 0.85 a:omg/omgth=1.12122，gam/gamth=1.61924
### test125 10 
slowing-down or bi-maxwell

### test126 10 9667794
q=1.75
B0=2T
R0=2.35m
a=0.47m
T_i=0.1keV
nrx=80
ntheta=16
nvpara=64
nmu=16
**output** 0.85 a:omg/omgth=1.14114，gam/gamth=1.42424

### test127 10 
carbon impurity only
$H^+,C^{6+}$
$A_{eff}=(n_iA_i+n_zA_z)/(n_i+n_z)A_i$
q=1.75
B0=2T
R0=2.35m
a=0.47m
R_R=0.9
L_R=0.8
T_i=0.1keV
nrx=80
ntheta=16
nvpara=64
nmu=16
Ionnum=2
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.8
  Upa    = 0.0
  Tem    = 0.10
/

&Par
  !Ion2
  mess   = 12.0
  charge = 6.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.033
  Upa    = 0.0
  Tem    = 0.10
  
| No. | n_i   | n_Z    | n_e   | job_id  | A_eff       |
|-----|-------|--------|-------|---------|-------------|
| 1   | 1.00  | 0.000  | 1.00  | 9667794 | 1           |
| 2   | 0.90  | 0.017  | 1.00  | 9668489 | 1.2         |
| 3   | 0.80  | 0.033  | 1.00  | 9668710 | 1.44        |
| 4   | 0.70  | 0.050  | 1.00  | 9671644 | 1.733333333 |
| 5   | 0.60  | 0.067  | 1.00  | 9671808 | 2.1         |
| 6   | 0.50  | 0.083  | 1.00  | 9671953 | 2.571428571 |
| 7   | 0.40  | 0.100  | 1.00  | 9672492 | 3.2         |
| 8   | 0.30  | 0.117  | 1.00  | 9672680 | 4.08        |
| 9   | 0.20  | 0.133  | 1.00  | 9673290 | 5.4         |
| 10  | 0.10  | 0.150  | 1.00  | 9673329 | 7.6         |
| 11  | 0.00  | 0.167  | 1.00  | 9673489 | 12          |

**output** 结果不符原因：环径比？
A_eff=12时PhiZ为NaN

### test128 10 9705575 digr=0.35
carbon impurity only
$H^+,C^{6+}$
$A_{eff}=(n_iA_i+n_zA_z)/(n_i+n_z)A_i$
q=1.75
B0=2T
R0=2.35m
a=0.47m
R_R=0.45
L_R=0.30
r_min=0.29
r_max=0.46
T_i=0.1keV
nrx=100
ntheta=16
nvpara=64
nmu=16
Ionnum=1
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.10
/

&Par
  !Ion2
  mess   = 12.0
  charge = 6.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.0
  Upa    = 0.0
  Tem    = 0.10

**output** 0.35 a:omg/omgth=1.15735，gam/gamth=1.18125
0.375 a:omg/omgth=1.16037，gam/gamth=1.05967
0.4 a:omg/omgth=1.15984，gam/gamth=1.13033

### test129 10 

^3688ec

carbon impurity only
$H^+,C^{6+}$
$A_{eff}=(n_iA_i+n_zA_z)/(n_i+n_z)A_i$
q=1.75
B0=2T
R0=2.35m
a=0.47m
R_R=0.45
L_R=0.30
r_min=0.29
r_max=0.46
T_i=0.1keV
nrx=$[0.17* 0.47/\rho_i]$   200
ntheta=16
nvpara=64
nmu=16
Ionnum=2
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.10
/

&Par
  !Ion2
  mess   = 12.0
  charge = 6.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.0
  Upa    = 0.0
  Tem    = 0.10

### test130 10 9708163
carbon impurity only
$H^+,C^{6+}$
$A_{eff}=(n_iA_i+n_zA_z)/(n_i+n_z)A_i$
q=3.5
B0=8.2T
R0=2.35m
a=0.47m
R_R=0.45
L_R=0.30
r_min=0.29
r_max=0.46
T_i=2keV
nrx=$[0.17* 0.47/\rho_i]$   110
ntheta=16
nvpara=64
nmu=16
Ionnum=1
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 2.0
/

&Par
  !Ion2
  mess   = 12.0
  charge = 6.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.0
  Upa    = 0.0
  Tem    = 0.10
**output** 与四阶FOW理论不符合

### 程序时间归一化 m_ref_new T_ref，需要修改timelist到主离子


### 高能粒子
#### double shifted Maxwellian distribution
$$F_{h0}=\frac{1}{2}\frac{n_h}{(2\pi T_h\big/m_h)^{3/2}}e^{-\frac{\mu B}{T_h}}\left(e^{-\frac{(v_\parallel-v_0)^2}{2T_h\big/m_h}}+e^{-\frac{(v_\parallel+v_0)^2}{2T_h\big/m_h}}\right)$$

[[#^fd13c9|EPs1.0]]
VlasovF0.f90 ComputInitF0
Global.f90 FDSHM
修改归一化因子m_ref_new,T_ref为主离子的参数，READparticle两遍
**EPs%PType=5**

V_0=3.5V_th -->（程序中归一化为C_sh sqrt(T_h/m_h)）4.9497C_sh
q=1.5
T_ep=T_i=T_e

#### Slowing-down distribution function

^976bb7

$$f_0=c\frac{n_0H(v_0-v)}{v^3+v_c^3}exp\left[-\left(\frac{\Lambda-\Lambda_0}{\Delta\Lambda}\right)^2\right]$$

$$c=\frac{3}{2\pi\left[\ln{(v_0^3+v_c^3)}-\ln v_c^3\right]\int_0^1\frac{1}{\sqrt{1-\Lambda}}exp\left[-\left(\frac{\Lambda-\Lambda_0}{\Delta\Lambda}\right)^2\right]\mathrm{d}\Lambda}$$
```Fortran
B = mesh_B(i,k)
    vtotal = sqrt((Parlist(iPar)%Vvpara0(l))**2+2*mu*B/ms)
    kinetener = 0.5_dp*ms*Parlist(iPar)%Vvpara0(l)**2+mu*B
    Lambda = mu*B0/kinetener
    if (U.ge.vtotal) then
      vc = sqrt(2.0_dp*14.8_dp*Parlist(0)%Tt(i)/ms)
      c = 2.0_dp*PI*(log(U**3+vc**3)-log(vc**3))*SDFfactor/3.0_dp
      c = 1.0_dp/c
      SDF = N*exp(-((Lambda-Lambda0sdf)/DeltaLambdasdf)**2)/(vtotal**3+(vc)**3)
      print*,"H1","B0:",B0,"B:",B
    else
      SDF = 0.0_dp
      print*,"H0","B0:",B0,"B:",B
    end if
```
v_0 birth speed -->Upa
PType = 6
critical energy[[zarzoso2012fully]]
$E_c=14.8A_hT_e\left(\frac{1}{n_e}\sum_j\frac{n_jZ_j^2}{A_j}\right)^{2/3}$

$\frac{1}{2}\int_0^{v_t^2-p}(\frac{-3}{p+c+\frac{t^3}{\sqrt{p+c}}}+\frac{c}{\frac{(p+c)^2}{4}}\times\frac{\frac{c}{c+p}-\lambda_0}{\delta\lambda^2})\frac{const*n_0*\sqrt{p}}{(c+p)^{3/2}+t^3}\times\exp{-(\frac{\frac{c}{p+c}-\lambda_0}{\delta\lambda})^2}\mathrm{d}c-\frac{const*n_0}{v_t^3+t^3}*\sqrt{p}*\exp{-(\frac{\frac{v_t^2-p}{v_t^2}-\lambda_0}{\delta\lambda})^2}$

equilibrium density distribution(芯部高，量级0.02* ne)
```fortran
DenN0 = exp(-kappaN*ra/R0*( r - rm &

                   - deltaN*tanh((r-Vmr(1)  )/deltaN) &

                   - deltaN*tanh((r-Vmr(nrx))/deltaN) ))
```
### test131 12 9730051
Upar=4
m_i=m_h=1
e=1
n_h=0.05
n_i=0.95
q=2
### test132 12
### test133 12 9735254

^45d464

Upar = 4
m_i=m_h=1
e=1
n_h=0.05
n_i=0.95
q=2
### test134 对比 12 [[#^45d464|test133]] 完全一致
#### try1 Upar=0  PType=5 9736423
#### try2 Upar=0 PType=4 9736448

### test135 10 Timelist  对比 
#### try1 m_ref_new=1 T_ref=0.1 9737424

#### try2 Timelist m_ref_new=2 T_ref=0.15 [[#^3688ec|test129 no.20]]

time1=time2. * sqrt(10* 2/15);
$T=t1\times\frac{R_0}{C_s^*}=t2\times\frac{R_0}{C_s}$
**结果不准确，m_ref_new和T_ref修改后Timelist反而更小**，时间增量和参考质量和温度有关，之前未使用主粒子归一化的在处理数据时均采用
time=time* sqrt(T* 2m_p /(m* 0.15keV))
处理

### test136 10 3165585
q=1.75
B0=2T
R0=2.35m
a=0.47m
R_R=0.55
L_R=0.45
r_min=0.44
r_max=0.56
T_i=0.1keV
T_e=0
nrx>$[0.12* 0.47/\rho_i]$   80
ntheta=16
nvpara=64
nmu=16
Ionnum=1
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.10
/
**output** 0.5 a:omg/omgth=1.00127，gam/gamth=1.0557
fit: omega=2.22328，gamma=-0.0716413
### bfun=fit(time(1:digtimeright),normE(1:digtimeright),fitf, 'StartPoint', initialGuess);%拟合效果

### Ntasks=16 和 mpirun -np 16 ./NLT >& ./run.out保持一致

有无I_MPI_PMI_LIBRARY定义 区别   3156674 无；   3156680 有 快
**相同参数两边跑的结果不一样？？** 大致一样，区别很小

### test137 10 3165593
q=1.75
B0=2T
R0=2.35m
a=0.47m
R_R=0.55
L_R=0.45
r_min=0.44
r_max=0.56
T_i=0.1keV
T_e=0
nrx>$[0.12* 0.47/\rho_i]$   80
ntheta=16
nvpara=64
nmu=16
Ionnum=2
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.88
  Upa    = 0.0
  Tem    = 0.10
/
&Par
  !Ion2
  mess   = 12.0
  charge = 6.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.02
  Upa    = 0.0
  Tem    = 0.10
**output** 0.5 a:omega=2.04253，gamma=-0.107321

### test138 12 [[#^45d464|test133]] 3166253
数据对比  bj th

### test139 10 3167651
q=1.75
B0=2T
R0=2.35m
a=0.47m
R_R=0.55
L_R=0.45
r_min=0.44
r_max=0.56
T_i=0.1keV
T_e=0   Tem    = 0.001
nrx>$[0.12* 0.47/\rho_i]$   80
ntheta=16
nvpara=64
nmu=16
Ionnum=2
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.76
  Upa    = 0.0
  Tem    = 0.10
/
&Par
  !Ion2
  mess   = 12.0
  charge = 6.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.04
  Upa    = 0.0
  Tem    = 0.10
**output** 0.5 a: omega=1.90182，gamma=-0.144516

### test140 10 3168843
q=1.75
B0=2T
R0=2.35m
a=0.47m
R_R=0.55
L_R=0.45
r_min=0.44
r_max=0.56
T_i=0.1keV
T_e=T_i=T_z
nrx>$[0.12* 0.47/\rho_i]$   80
ntheta=16
nvpara=64
nmu=16
Ionnum=2
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.76
  Upa    = 0.0
  Tem    = 0.10
/
&Par
  !Ion2
  mess   = 12.0
  charge = 6.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.04
  Upa    = 0.0
  Tem    = 0.10
  **output** 0.5 a: omega=2.26878，gamma=-0.0462198

### test141 10 3168856
q=1.75
B0=2T
R0=2.35m
a=0.47m
R_R=0.55
L_R=0.45
r_min=0.44
r_max=0.56
T_i=0.1keV
T_e=T_i=T_z
nrx>$[0.12* 0.47/\rho_i]$   80
ntheta=16
nvpara=64
nmu=16
Ionnum=2
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.88
  Upa    = 0.0
  Tem    = 0.10
/
&Par
  !Ion2
  mess   = 12.0
  charge = 6.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.02
  Upa    = 0.0
  Tem    = 0.10
  **output** 0.5 a: omega=2.38926，gamma=-0.0368397

### test142 10 3168867
q=1.75
B0=2T
R0=2.35m
a=0.47m
R_R=0.55
L_R=0.45
r_min=0.44
r_max=0.56
T_i=0.1keV
T_e=T_i=T_z
nrx>$[0.12* 0.47/\rho_i]$   80
ntheta=16
nvpara=64
nmu=16
Ionnum=1
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.10
/
**output** 0.5 a: omega=2.56254，gamma=-0.0540856
### test143 10 3169348
q=1.75
B0=2T
R0=4m
a=0.47m
R_R=0.55
L_R=0.45
r_min=0.44
r_max=0.56
T_i=0.1keV
T_e=T_i=T_z
nrx>$[0.12* 0.47/\rho_i]$   80
ntheta=16
nvpara=64
nmu=16
Ionnum=1
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.10
/
**output**   0.5 a: omega=2.57403，gamma=-0.0550123

### test144 10 3169358
q=1.75
B0=2T
R0=4m
a=0.47m
R_R=0.55
L_R=0.45
r_min=0.44
r_max=0.56
T_i=0.1keV
T_e=T_i=T_z
nrx>$[0.12* 0.47/\rho_i]$   80
ntheta=16
nvpara=64
nmu=16
Ionnum=2
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.88
  Upa    = 0.0
  Tem    = 0.10
/
&Par
  !Ion2
  mess   = 12.0
  charge = 6.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.02
  Upa    = 0.0
  Tem    = 0.10
**output**  0.5 a: omega=2.39762，gamma=-0.0382555

### test145 10 3169370
q=1.75
B0=2T
R0=4m
a=0.47m
R_R=0.55
L_R=0.45
r_min=0.44
r_max=0.56
T_i=0.1keV
T_e=T_i=T_z
nrx>$[0.12* 0.47/\rho_i]$   80
ntheta=16
nvpara=64
nmu=16
Ionnum=2
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.76
  Upa    = 0.0
  Tem    = 0.10
/
&Par
  !Ion2
  mess   = 12.0
  charge = 6.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.04
  Upa    = 0.0
  Tem    = 0.10
**output**  0.5 a: omega=2.27379，gamma=-0.0433433

### test146 12 [[#^45d464|test133]] 3172228
q=2
B0=2T
R0=2.35m
a=0.47m
R_R=0.55
L_R=0.45
r_min=0.44
r_max=0.56
T_i=0.1keV
nrx=80
ionnum=2
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.10
/

&Par
  !Ion2
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.0
  Tem    = 0.10
/


### test147 12 1079353
q=2
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.1keV
nrx=120
ionnum=2
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.10
/

&Par
  !Ion2
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.0
  Tem    = 0.10
/

### test148 12 3190242
q=2
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.1keV
nrx=120
ionnum=2
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.97
  Upa    = 0.0
  Tem    = 0.10
/

&Par
  !Ion2
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.03
  Upa    = 4.0
  Tem    = 0.10
/

### test149 12 3190828
q=2
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.1keV
nrx=120
ionnum=2
&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.10
/

&Par
  !Ion2
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 2.83
  Tem    = 0.10
/

### test150 12 3191353
q=2
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.1keV
nrx=120
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.10
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.00
  Tem    = 0.10
/

### test151 12 3194715
q=2
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.1keV
nrx=120
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.10
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.00
  Tem    = 0.20
/
**output** 出现不稳定的模，可能为EGAM中GAM演变而来的分支  0.525 a: omega=1.98448，gamma=0.160148

### test152 12 3202430
q=2
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.08keV
nrx=100
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.08
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.00
  Tem    = 0.20
/
**output** 
0.525 a: omega=2.04958，gamma=0.176025
### test153 12 1092733
q=2
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.00
  Tem    = 0.20
/
**output** 0.525 a: omega=2.20512，gamma=0.0447751
### test154 12 1093300  [reference]
q=1
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 2.00
  Tem    = 0.06
/


### test155 12 1094610
q=1
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=200
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 2.00
  Tem    = 0.06
/

### test156 12 1094672
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 2.00
  Tem    = 0.06
/

### test157 12 1094703
q=1
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 3.00
  Tem    = 0.06
/

### test158 12 1101143
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 3.00
  Tem    = 0.06
/
**output** 出现两支阻尼的模
#### 158-1time* 2  1106682

### test159 12 1101787
q=2.0
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.9497
  Tem    = 0.03
/

### test160 12 1102136 符合EGAM理论
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.9497
  Tem    = 0.03
/
重复1106643
**output** 0.525 a:omg/omgth=1.04387，gam/gamth=0.882578
fit: omega=2.44864，gamma=0.201143
### test161 12 1107910

^dba2f7

q=2.7
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.2426
  Tem    = 0.03
/
**output** 0.525 a: omega=2.1444，gamma=-0.011285
与理论不符合
### test162 12 1110702

^3fc273

q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.2426
  Tem    = 0.03
/
**output** 0.525 a:omg/omgth=1.08221，gam/gamth=0.847791
fit: omega=2.3899，gamma=0.118161
### test163 12 1116022 [[#^3fc273|test162]]
time* 2
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.2426
  Tem    = 0.03
/

### test164 12 1128252
q=3.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.2426
  Tem    = 0.03
/
**output** 相比于[[#^dba2f7|test161]]，q增加却由阻尼变为增长
### test165 12 1160557
q=4.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.2426
  Tem    = 0.03
/
**output** 0.525 a: omega=2.06438，gamma=0.000472468
### test166 12 1161275
q=4.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 3
  Tem    = 0.03
/
**output** 0.525 a: omega=2.39031，gamma=-0.00722475
theory 2.75056 -0.0387327 i
### test167 12 1165168
q=4.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 2.5
  Tem    = 0.03
/

### test168 12 1207173
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  [vMax   = 10.0]
  PType  = 5
  Den    = 0.10
  Upa    = 4.9497
  Tem    = 0.03
/

### 改变vMax：看代码 F0--vpara

### test169 SDFtest 12399-->1232593
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 6
  Den    = 0.10
  Upa    = 3.60
  Tem    = 0.03
/
Lambda0sdf = 0.50
DeltaLambdasdf = 0.05
```

### test170 SDFtest 1233906

^b93f55

q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 6
  Den    = 0.10
  Upa    = 3.60
  Tem    = 0.03
/
Lambda0sdf = 0.50
DeltaLambdasdf = 0.10
```

### test171 SDFtest 1233921
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 6
  Den    = 0.10
  Upa    = 3.60
  Tem    = 0.03
/
Lambda0sdf = 0.50
DeltaLambdasdf = 0.15
```

### test172 SDFtest 1233998
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 6
  Den    = 0.10
  Upa    = 3.60
  Tem    = 0.03
/
Lambda0sdf = 0.50
DeltaLambdasdf = 0.20
```

### test173 SDFtest 1234076
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 6
  Den    = 0.10
  Upa    = 3.60
  Tem    = 0.03
/
Lambda0sdf = 0.50
DeltaLambdasdf = 0.02
```

### test174 10 12 SDFtest GAM结果对比
```fortran
r_max = 0.61
r_min = 0.44
nrx = 150
R_R = 0.60
L_R = 0.45
Ionnum=1
T_e = 0.03
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.00
  Upa    = 0.0
  Tem    = 0.03
/
q=1.5
B0=2T
R0=2.35m
a=0.47m
10-1234392 12-1234448 SDFtest-1234477
```
**output** theory:2.76402 -0.0317919 i
simulation 2.72201-0.00901i
### test175 SDFtest 1235448
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6  
  Den    = 0.10
  Upa    = 4.86 
  Tem    = 0.03
/
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
```
**output**：两支模，低频支增长，高频支阻尼
0.525 a: omega1=1.75378，gamma1=0.0294815;
omega2=3.11168，gamma2=-0.00861469
### test176 SDFtest 1239904
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.08
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6  
  Den    = 0.10
  Upa    = 4.86 
  Tem    = 0.08
/
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.08keV
T_e=0.08kev
nrx=100
ionnum=2
```
**output**:0.525 a: omega1=3.0964，gamma1=0.0160669;omega2=1.7249，gamma2=-0.024318
### test177 10 1239954
```fortran
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.08keV
T_e=0.08kev
nrx=100
ionnum=1
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.08
/       
```
**output**:theory2.75461 -0.0436017i
0.525 a: omega=2.70765，gamma=-0.021786
### test178 SDFtest 1250819
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.80
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6
  Den    = 0.20
  Upa    = 9.72
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.05
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
```

### test179 SDFtest 1251126
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.80
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6
  Den    = 0.20
  Upa    = 9.72
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.01
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
```

### test180 SDFtest 1251368

^a3313f

```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.01
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
```
**类似两支不稳定的模式**

### 输出归一化 C_s/R -> V_ti/R
### test181 SDFtest 1269945
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6
  Den    = 0.10
  Upa    = 9.72
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.01
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
```
**output**: 0.525 a: omega=1.61389，gamma=0.0759483
### test182 13 [[#test180 SDFtest 1251368|test80]] 1270185
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.01
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
FlagKE = .F. -> .T.
dtMax = 2.0-> 0.1  !maxstep=100000**停止**
```

### test183 13 1271498
```Fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.01
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
FlagKE = .F. -> .T.
dtMax = 2.0-> 0.1  
maxstep=200000 !100000->200000
```

### test184 13 1283025 [[#^a3313f|test180]]
输出分布函数网格，大小
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.01
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
FlagKE = .F.
dtMax = 2.0
maxstep=100000 !100000
nvpara = 64
        muType = 3
           nmu = 16
```

### test185 13 1285223

^ca0579

```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.01
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
FlagKE = .F.
dtMax = 2.0
maxstep=100000 !100000
nvpara = 128
        muType = 3
           nmu = 16
```
**nvpara** 64->128  一支增长->两支增长?(qk_r rho_i 0.1496) [[#test180 SDFtest 1251368|test180]]
### test186 13 1285887
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 3.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.01
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
FlagKE = .F.
dtMax = 2.0
maxstep=100000 !100000
nvpara = 128
        muType = 3
           nmu = 16
```

### t187 13 1285897
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.05
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
FlagKE = .F.
dtMax = 2.0
maxstep=100000 !100000
nvpara = 128
        muType = 3
           nmu = 16
```

### pro_id-->0--nmu-1

### t188 13 [[#^ca0579|test185]] 1286079
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 3.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.01  !delta pitch angle 稍微增加
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
FlagKE = .F.
dtMax = 2.0
maxstep=100000 !100000
nvpara = 96
        muType = 3
           nmu = 16
```
Vpara从4改为3，网格密度不变，振荡模式区别大，4覆盖的范围广
### t189 13  1288515
不同进程的F0不相同，即mu_loc
输出数据应该是nmu* nvpara行数据
```fortran
pro_id=0-15!输出速度网格、F0  ，不同进程的F0相同，即mu_loc
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.05
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
FlagKE = .F.
dtMax = 2.0
maxstep=100000 !100000
nvpara = 128
        muType = 3
           nmu = 16

```

### t190 13 1289105
```fortran
pro_id=0-15!输出速度网格、F0  ，不同进程的F0相同，即mu_loc
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.05
R0 = 2.35
         a = 0.47
      nofr = 513                    !radial mesh points
  noftheta = 257                    !poloidal mesh points
     Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
        e0 = 0.180776594471700d0
        q0 = 1.4 !1.455125d0
        s0 = 0.8 !0.826217678893566d0
      F_cn = T
        c0 =  1.50
    
Ionnum=2
FlagKE = .F.
dtMax = 2.0
maxstep=100000 !100000
nvpara = 128
        muType = 3
           nmu = 32
```
**nmu** 32和16震荡模式区别不大 

### t191 13  [[#test185 13 1285223|test185]] 1290493
**vMax**
pro_id 输出网格、F0
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 5.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.01
R0 = 2.35
a = 0.47
Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
c0 =  1.50
Ionnum=2
FlagKE = .F.
dtMax = 2.0
maxstep=100000 !100000
nvpara = 160
muType = 3
nmu = 16
```

### t192 13 [[#t191 13 test185 13 1285223 test185 1290493|t191]] 1290981
**vMax**
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 6.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.01
R0 = 2.35
a = 0.47
Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
c0 =  1.50
Ionnum=2
FlagKE = .F.
dtMax = 2.0
maxstep=100000 !100000
nvpara = 192
muType = 3
nmu = 16
```

### t193 13 [[#t192 13 t191 13 test185 13 1285223 test185 1290493 t191 1290981|t192]] 1291259
**vMax**
```fortran
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/       

&Par         
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 7.0 
  PType  = 6
  Den    = 0.10
  Upa    = 4.86
  Tem    = 0.12
/
&Par
  !Electron
  Tem    = 0.12
/ 
Lambda0sdf = 0.50
    DeltaLambdasdf = 0.01
R0 = 2.35
a = 0.47
Bphi0 = 2.0                    !toroidal magnetic at axis (Tesla)
c0 =  1.50
Ionnum=2
FlagKE = .F.
dtMax = 2.0
maxstep=100000 !100000
nvpara = 224
muType = 3
nmu = 16
```

mu=(1:nmu)0.4981 vMax^2    m/2B=0.**4981**
### t194 13 1291865
vmax=5,nmu=16* 3
nvpara=160
DeltaLambdasdf = 0.01

### t195 13 1293832
vmax=5,nmu=16* 3
nvpara=160
DeltaLambdasdf = 0.05

### t196 13 1293846
vmax=5,nmu=16* 2
nvpara=160
DeltaLambdasdf = 0.05

### t197 13 1294720
vmax=5,nmu=16* 1
nvpara=160
DeltaLambdasdf = 0.05

### t198 13 1294859 [[#test160 12 1102136 符合EGAM理论|test160 PT=5]]
```fortran
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.9497
  Tem    = 0.03
/
```

### t199 13 1294908 [[#t197 13 1294720|t197]]
vmax=4,nmu=16* 2
nvpara=128
DeltaLambdasdf = 0.05

### t200 13 1295146 [[#test160 12 1102136 符合EGAM理论|test160 PT=5]]
```fortran
q=1.5
B0=2T
R0=2.35m
a=0.47m
R_R=0.60
L_R=0.45
r_min=0.44
r_max=0.61
T_i=0.03keV
T_e=0.03kev
nrx=150
ionnum=2
&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.03
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 9.0
  PType  = 5
  Den    = 0.10
  Upa    = 4.9497
  Tem    = 0.03
/
nvpara=144
nmu=96

```


## test201-
### t201 13 [[#t196 13 1293846|t196]] 1295318
nmu初步确定为32，vMax =5
```fortran
vmax=5,nmu=16* 2
nvpara=240
DeltaLambdasdf = 0.05
```

### t202 13 1295711
```fortran
vmax=5,nmu=16* 2
nvpara=80
DeltaLambdasdf = 0.05
```

### t203 13 1295889
```fortran
vmax=5,nmu=16
nvpara=80
DeltaLambdasdf = 0.05
```

### slowing down velocity grid
```Fortran
vmax=5,nmu=16* 2
nvpara=160
DeltaLambdasdf = 0.05
Upa = 4.86
```

### t204 13 
scan Upa
```Fortran
vmax=5,nmu=16* 2
nvpara=160
DeltaLambdasdf = 0.05
Upa = 3:0.5:6
Te=Ti=Th=120ev
nh=0.1
m=2
e=1
q=1.5
B=2
```
拟合经验：拟合区域调为两波叠加区域
V0大小，分布改变  关系:V_0=6,vMax=5√

| V_0(Upa) | job_id          | omega/ vti/R | gamma/ vti/R | （omega2）          | （gamma2）           |
|----------|-----------------|--------------|--------------|-----------|------------|
| 3        | 1302280.000000  | 1.831410     | -0.033211    |           |            |
| 3.5      | 1302283.000000  | 1.907380     | 0.013300     |           |            |
| 4        | 1302286.000000  | 1.990870     | 0.020039     |           |            |
| 4.5      | 1302289.000000  | 2.087750     | 0.010676     | 1.179990  | 0.008626   |
| 5        | 1302292.000000  | 2.197250     | 0.007331     | 1.245390  | -0.007002  |
| 5.5      | 1302294.000000  | 2.315200     | 0.005178     | 1.335990  | 0.039297   |
| 6        | 1302296.000000  | 2.441930     | 0.014462     | 1.391820  | 0.057781   |

### t205 13 1303534
```fortran
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.12
/

&Par
  !Ion1
  mess   = 2.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.90
  Upa    = 0.0
  Tem    = 0.12
/

&Par
  !Ion2
  mess   = 2.0
  charge = 1.0
  vMax   = 9.0
  **PType  = 5**
  Den    = 0.10
  Upa    = 4.5
  Tem    = 0.12
/
nrx=80
ntheta=16
nvpara=144
nmu=96
B=2
q=1.5
```
### debug slowing down EP
critienerfasdf计算错误，**Fortran 整数除法中加上. _ dp 幂次**
E_0 !=T_h
### t206 13 1327177 ##二阶FOW共振证据
```fortran
vmax=5,nmu=16* 2
nvpara=160
DeltaLambdasdf = 0.05
Upa = 5.5
Lambda0=0.5
Te=Ti=Th=120ev
nh=0.1
m=2
e=1
q=1.5
B=2
```


### t207 13 1330312
```fortran
R0=7.05
a=0.47
vmax=5,nmu=16* 2
nvpara=160
DeltaLambdasdf = 0.05
Upa = 5.5
Lambda0=0.5
Te=Ti=Th=120ev
nh=0.1
m=2
e=1
q=1.5
B=2
```
**output**
(2.40245 -2.33303)/2.40245 
0.0289
 (0.00559796-0.00480436)/0.00559796
0.1418
(1.3414640670247546-1.3231)/1.3414640670247546
0.0137
 (0.07416155516435373-0.0575742)/0.07416155516435373
0.2237

### t208 13 1330316
```fortran
R0=7.05
a=0.47
vmax=5,nmu=16* 2
nvpara=160
DeltaLambdasdf = 0.02
Upa = 5.5
Lambda0=0.5
Te=Ti=Th=120ev
nh=0.1
m=2
e=1
q=1.5
B=2
```
(2.40245-2.33686)/2.40245 
0.0273
(0.00802091-0.00559796)/0.00559796
0.4328
(1.3414640670247546-1.33191)/1.3414640670247546
0.0071
(0.07416155516435373-0.0636559)/0.07416155516435373
0.1417


### t209 13 1335463
```fortran
R0=7.05
a=1.5
m_h=2
m_i=2
q=2
Lambda0=0.3
R_L=0.10
R_R=0.80
rmin=0.09
rmax=0.81
T_e=4kev
T_i=4kev
T_h=340 000ev
n_h=0.25
n_i=0.75
nrx=180
deltaLambda=0.05
e=1
Upa=1
B=2
vMax=4
```

### t210 13 1335585
```fortran
R0=7.05
a=1.5
m_h=2
m_i=2
q=2
Lambda0=0.3
R_L=0.10
R_R=0.80
rmin=0.09
rmax=0.81
T_e=4kev
T_i=4kev
T_h=340 000ev
n_h=0.25
n_i=0.75
nrx=180
deltaLambda=0.05
e=1
Upa=1
B=2
vMax=0.4
```

### t211 13 1335690
```fortran
R0=7.05
a=1.5
m_h=2
m_i=2
q=2
Lambda0=0.3
R_L=0.10
R_R=0.80
rmin=0.09
rmax=0.81
T_e=4kev
T_i=4kev
T_h=340 000ev
n_h=0.25
n_i=0.75
nrx=180
deltaLambda=0.05
e=1
Upa=1
B=2
vMax=0.6
```

### t212 13 1340642
```fortran
R0=7.05
a=1.5
m_h=2
m_i=2
q=2
Lambda0=0.3
R_L=0.10
R_R=0.80
rmin=0.09
rmax=0.81
T_e=4kev
T_i=4kev
T_h=340 000ev
n_h=0.25
n_i=0.75
nrx=180
deltaLambda=0.05
e=1
Upa=1
B=2
vMax=0.9
```

### t213 13 1381402
```fortran
R0=7.05
a=1.5
m_h=2
m_i=2
q=2
Lambda0=0.3
R_L=0.10
R_R=0.80
rmin=0.09
rmax=0.81
T_e=4kev
T_i=4kev
T_h=340 000ev
n_h=0.25
n_i=0.75
nrx=260
deltaLambda=0.05
e=1
Upa=1
B=3
vMax=0.9
```

### t214 13 1381409
```fortran
R0=7.05
a=0.47
m_h=2
m_i=2
q=2
Lambda0=0.3
R_L=0.10
R_R=0.80
rmin=0.09
rmax=0.81
T_e=4kev
T_i=4kev
T_h=340 000ev
n_h=0.25
n_i=0.75
nrx=80
deltaLambda=0.05
e=1
Upa=1
B=3
vMax=0.9
```
**output** 两个频率相差 0.3V_ti/R的波

### t215 13 1381711
```fortran
R0=7.05
a=0.47
m_h=2
m_i=2
q=2
Lambda0=0.3
R_L=0.10
R_R=0.80
rmin=0.09
rmax=0.81
T_e=4kev
T_i=4kev
T_h=340 000ev
n_h=0.25
n_i=0.75
nrx=80
deltaLambda=0.05
e=1
Upa=1
B=3
vMax=0.9
simTotalTime = 80.4
	       maxStep = 200000
```


### t216 13 1383615
```fortran
R0=7.05
a=0.47
m_h=2
m_i=2
q=2
Lambda0=0.3
R_L=0.10
R_R=0.80
rmin=0.09
rmax=0.81
T_e=4kev
T_i=4kev
T_h=340 000ev
n_h=0.25
n_i=0.75
nrx=80
deltaLambda=0.05
e=1
Upa=1
B=3
vMax=0.9
simTotalTime = 120.6
       maxStep = 300000
nvpara=128
nmu=48
FlagKE=.F.
dtMax=2.0
```

### t217 13 1413723
```fortran
FlagKE=.T.
dtMax=0.1
R0=7.05
a=0.47
m_h=2
m_i=2
q=2
Lambda0=0.3
R_L=0.10
R_R=0.80
rmin=0.09
rmax=0.81
T_e=4kev
T_i=4kev
T_h=340 000ev
n_h=0.25
n_i=0.75
nrx=80
deltaLambda=0.05
e=1
Upa=1
B=3
vMax=0.9
simTotalTime = 120.6
       maxStep = 300000
nvpara=128
nmu=48
```

### LHD experimental condition[[ido2015identification|article]]
R0=3.75m
a=0.6m
B0=1.375T
n0=$0.1\times10^{19}m^{-3}$

### t218 13 1421529
```fortran
n_h=0.06
T_i=0.6 KeV
T_e=4 KeV
E0=175 KeV
\Lambda_0=0.116
Delta-Lambda=0.1
q=2.5
Delta-L=0.3
inside 0.6a
R0=3.75m
a=0.6m
B0=1.375T
n0=$1\times10^{19}m^{-3}$
m_i=m_h=2
e_i=e_h=1
0.2--0.5
vmax=0.9  fast ion
```

| 参量（理论对比） |  |  |  |
| ---- | ---- | ---- | ---- |
| E_0/T_e | 43.7500 | lambda_i | 0.1271 |
| E_c/T_e | 18.6468 | lambda_h | 2.0404 |
|  |  |  |  |
| tau_c | 0.4262 |  |  |
|  |  | nrx | 52.7418 |
| E_c/KeV | 74.5873 |  |  |
|  |  |  |  |
| n_i/n_e | 0.9400 |  |  |
| 参量（程序调参） |  |  |  |
| T_h/KeV | 350 |  |  |
| v_0/C_s^h | 1 |  |  |

### t219 13 1421531
```fortran
n_h=0.06
T_i=0.6 KeV
T_e=4 KeV
E0=175 KeV
\Lambda_0=0.116
Delta-Lambda=0.1
q=2.5
Delta-L=0.3
inside 0.6a
R0=3.75m
a=0.6m
B0=1.375T
n0=$0.1\times10^{19}m^{-3}$
m_i=m_h=2
e_i=e_h=1
0.2--0.5
vmax=0.9  fast ion
```

| 参量（理论对比） |  |  |  |
| ---- | ---- | ---- | ---- |
| E_0/T_e | 43.7500 | lambda_i | 0.1271 |
| E_c/T_e | 18.6468 | lambda_h | 2.0404 |
|  |  |  |  |
| tau_c | 0.4262 |  |  |
|  |  | nrx | 52.7418 |
| E_c/KeV | 74.5873 |  |  |
|  |  |  |  |
| n_i/n_e | 0.9400 |  |  |
| 参量（程序调参） |  |  |  |
| T_h/KeV | 350 |  |  |
| v_0/C_s^h | 1 |  |  |
### t220 13 22292
```fortran
n_h=0.06
T_i=0.6 KeV
T_e=4 KeV
E0=175 KeV
\Lambda_0=0.116
Delta-Lambda=0.1
q=2.5
Delta-L=0.3
inside 0.6a
R0=3.75m
a=0.6m
B0=1.375T
n0=$1\times10^{19}m^{-3}$
m_i=m_h=2
e_i=e_h=1
0.2--0.5
vmax=1.5  fast ion
```

### t221 13 1424063
```fortran
n_h=0.06
T_i=0.6 KeV
T_e=4 KeV
E0=175 KeV
\Lambda_0=0.116
Delta-Lambda=0.1
q=2.5
Delta-L=0.3
inside 0.6a
R0=3.75m
a=0.6m
B0=1.375T
n0=$1\times10^{19}m^{-3}$
m_i=m_h=2
e_i=e_h=1
0.2--0.5
vmax=1.5  fast ion
Flag_KE T
dtmax =0.1
maxstep=200000
```
**output** 动理学电子 阻尼作用
### t222 13 1424918
```fortran
n_h=0.06
T_i=0.6 KeV
T_e=4 KeV
E0=175 KeV
\Lambda_0=0.116
Delta-Lambda=0.1
q=2.5
Delta-L=0.5
inside 0.6a
R0=3.75m
a=0.6m
B0=6.5 T
n0=$1\times10^{19}m^{-3}$
m_i=m_h=2
e_i=e_h=1
0.1--0.6
vmax=1.5  fast ion
Flag_KE F
dtmax =2.0
maxstep=100000
nrx=410
nvpara = 160
nmu = 32
```

| 参量（理论对比） |  |  |  |
| ---- | ---- | ---- | ---- |
| E_0/T_e | 43.7500 | lambda_i | 0.0161 |
| E_c/T_e | 18.6468 | lambda_h | 0.2590 |
|  |  |  |  |
| tau_c | 0.4262 |  |  |
|  |  | nrx | 405.1526 |
| E_c/KeV | 74.5873 |  |  |
|  |  |  |  |
| n_i/n_e | 0.9400 |  |  |
| 参量（程序调参） |  |  |  |
| T_h/KeV | 350 |  |  |
| v_0/C_s^h | 1 |  |  |
|  |  |  |  |

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|-0.00363|7.40622|0.21436|2.59881|
|模拟|-0.05602|7.20996|0.234103|2.5533|
|误差|1444.47%|2.65%|9.21%|1.75%|
### dtmax2.0-->0.1 time* 20


### t224 13 

^354b34

E_0=150kev  1441673

E_0=100kev  1441727



#### E_0=75kev  1441745

^7b6efb

DeltaLambdasdf = 0.1
R_R = 0.60 
L_R = 0.10 
ntheta = 16
nvpara = 160
nmu = 32

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 2.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 2.00 |  |  |  | E_0/T_e | 18.7500 | lambda_i | 0.0416 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 18.6468 | lambda_h | 0.1695 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.50 |  | tau_c | 0.9945 |  |  |
|  |  |  |  |  |  |  | nrx | 156.9149 |
| E_0/KeV | 75.00 | e_i | 1.00 |  | E_c/KeV | 74.5873 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 6.5000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.50 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 3750 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.002327 | 2.20563 | 0.074499 | 1.23638 |
| 模拟 | 0.014261 | 2.18276 | 0.070128 | 1.20534 |
| 误差 | 512.82% | 1.04% | 5.87% | 2.51% |

### t225 13 
**随delta lambda 增加，增长率增大**

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 2.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 2.00 |  |  |  | E_0/T_e | 18.7500 | lambda_i | 0.0416 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 18.6468 | lambda_h | 0.1695 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.50 |  | tau_c | 0.9945 |  |  |
|  |  |  |  |  |  |  | nrx | 156.9149 |
| E_0/KeV | 75.00 | e_i | 1.00 |  | E_c/KeV | 74.5873 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 6.5000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.50 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 3750 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |

[[#^7b6efb|Δλ=0.1]]  

Δλ=0.075：

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|0.002327|2.20563|0.074499|1.23638|
|模拟|0.013407|2.18895|0.0701|1.20961|
|误差|476.16%|0.76%|5.90%|2.17%|

Δλ=0.050：

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|0.002327|2.20563|0.074499|1.23638|
|模拟|0.009851|2.19505|0.070024|1.21485|
|误差|323.35%|0.48%|6.01%|1.74%|
Δλ=0.025 ×
Δλ=0.010 ×

Δλ=0.040:

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|0.002327|2.20563|0.074499|1.23638|
|模拟|0.003402|2.19416|0.066345|1.2295|
|误差|46.21%|0.52%|10.94%|0.56%|

Δλ=0.045:

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|0.002327|2.20563|0.074499|1.23638|
|模拟|0.007291|2.19553|0.068831|1.21942|
|误差|213.32%|0.46%|7.61%|1.37%|
Δλ=0.041:
Δλ=0.039:
Δλ=0.038:
Δλ=0.037:
**磁场参数不正确**

Δλ=0.035:
拟合不准确


Δλ=0.030:拟合不准确


### t226 13 1446611
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 18.7500 | lambda_i | 0.0294 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.1199 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.50 |  | tau_c | 0.7893 |  |  |
|  |  |  |  |  |  |  | nrx | 221.9112 |
| E_0/KeV | 75.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 6.5000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.50 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 3750 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
理论解
 2.20409 - 0.00165517
2.2041 + 0.00101744
1.23704 + 0.0639169

模拟结果 Δλ=0.05
2.22541  -0.00407537
1.21541+0.0606462

模拟结果 Δλ=0.0？

### t227 13 1446630
GAM 无高能离子

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 18.7500 | lambda_i | 0.0294 |
|  |  | n_h/n_e | 0.00 |  | E_c/T_e | 14.8000 | lambda_h | 0.1199 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.50 |  | tau_c | 0.7893 |  |  |
|  |  |  |  |  |  |  | nrx | 221.9112 |
| E_0/KeV | 75.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 1.0000 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 6.5000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.50 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 3750 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
fit: omega=1.73561，gamma=-0.0059094
theory 2thFOW FLR   1.75862 -0.00320748
1.31%                        84.24%

### t228 13 
LHD parameter
B=2T
GAM:3428160
EGAM:3428151

### t229 13 3429607
delta_lambda=0.045

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.8000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 18.7500 | lambda_i | 0.0957 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.1853 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.50 |  | tau_c | 0.7893 |  |  |
|  |  |  |  |  |  |  | nrx | 68.2804 |
| E_0/KeV | 75.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 2.0000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.50 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 3750 |  |  |
| 常量 |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
| mp | 1.6726E-27 |  |  |  |  |  |  |  |
| e0 | 1.6022E-19 |  |  |  |  |  |  |  |

### t230 13 3429849
delta lambda=0.05
0.04-0.86

| 输入参数                     |            |          |         |  |           |          |          |           |
|--------------------------|------------|----------|---------|--|-----------|----------|----------|-----------|
| m_h                      | 2.00       | Lambda 0 | 0.1160  |  | 参量（理论对比）  |          |          |           |
| m_i                      | 2.00       |          |         |  | E_0/T_e   | 43.7500  | lambda_i | 0.0477    |
|                          |            | n_h/n_e  | 0.06    |  | E_c/T_e   | 18.6468  | lambda_h | 0.7652    |
| T_e/KeV                  | 4.00       |          |         |  |           |          |          |           |
| T_i/KeV                  | 0.60       | q        | 2.50    |  | tau_c     | 0.4262   |          |           |
|                          |            |          |         |  |           |          | nrx      | 135.1508  |
| E_0/KeV                  | 175.00     | e_i      | 1.00    |  | E_c/KeV   | 74.5873  |          |           |
|                          |            | e_h      | 1.00    |  |           |          |          |           |
| a                        | 0.60       |          |         |  | n_i/n_e   | 0.9400   |          |           |
| R                        | 3.75       |          |         |  |           |          |          |           |
| B0/T                     | 1.3750     |          |         |  |           |          |          |           |
|                          |            |          |         |  |           |          |          |           |
| perturbation wave length | 0.80       |          |         |  | 参量（程序调参）  |          |          |           |
|                          |            |          |         |  | T_h/KeV   | 8750     |          |           |
| 常量                       |            |          |         |  | v_0/C_s^h | 0.2      |          |           |
| mp                       | 1.6726E-27 |          |         |  |           |          |          |           |
| e0                       | 1.6022E-19 |


### t231 13 3429851
delta lambda=0.05
0.04-0.86

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 43.7500 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.5411 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 0.3383 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 175.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8750 |  |  |
| 常量 |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
| mp | 1.6726E-27 |  |  |  |  |  |  |  |
| e0 | 1.6022E-19 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
#### 3429996

### t232 13 3429997
delta lambda=0.05
0.04-0.86

| 输入参数                     |            |          |         |  |           |          |          |           |
|--------------------------|------------|----------|---------|--|-----------|----------|----------|-----------|
| m_h                      | 1.00       | Lambda 0 | 0.1160  |  | 参量（理论对比）  |          |          |           |
| m_i                      | 1.00       |          |         |  | E_0/T_e   | 25.0000  | lambda_i | 0.0337    |
|                          |            | n_h/n_e  | 0.06    |  | E_c/T_e   | 14.8000  | lambda_h | 0.4090    |
| T_e/KeV                  | 4.00       |          |         |  |           |          |          |           |
| T_i/KeV                  | 0.60       | q        | 2.50    |  | tau_c     | 0.5920   |          |           |
|                          |            |          |         |  |           |          | nrx      | 191.1320  |
| E_0/KeV                  | 100.00     | e_i      | 1.00    |  | E_c/KeV   | 59.2000  |          |           |
|                          |            | e_h      | 1.00    |  |           |          |          |           |
| a                        | 0.60       |          |         |  | n_i/n_e   | 0.9400   |          |           |
| R                        | 3.75       |          |         |  |           |          |          |           |
| B0/T                     | 1.3750     |          |         |  |           |          |          |           |
|                          |            |          |         |  |           |          |          |           |
| perturbation wave length | 0.80       |          |         |  | 参量（程序调参）  |          |          |           |
|                          |            |          |         |  | T_h/KeV   | 5000     |          |           |
| 常量                       |            |          |         |  | v_0/C_s^h | 0.2      |          |           |
| mp                       | 1.6726E-27 |          |         |  |           |          |          |           |
| e0                       | 1.6022E-19 |

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|-0.0092|5.59986|0.241652|2.52237|
|模拟|-0.01443|5.44638|0.198992|2.507|
|误差|56.85%|2.74%|17.65%|0.61%|
### t233 13 3430576
delta lambda=0.05
0.04-0.86

| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.1160 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 18.7500 | lambda_i | 0.0337   |
|                          |        | n_h/n_e  | 0.06   |     | E_c/T_e   | 14.8000 | lambda_h | 0.3542   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 0.60   | q        | 2.50   |     | tau_c     | 0.7893  |          |          |
|                          |        |          |        |     |           |         | nrx      | 191.1320 |
| E_0/KeV                  | 75.00  | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 0.60   |          |        |     | n_i/n_e   | 0.9400  |          |          |
| R                        | 3.75   |          |        |     |           |         |          |          |
| B0/T                     | 1.3750 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.80   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 3750    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |
**LHD,高频阻尼，低频增长**

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|-0.01014|5.02686|0.253776|2.42982|
|模拟|-0.01276|4.911|0.202478|2.39017|
|误差|25.87%|2.30%|20.21%|1.63%|
### t234 13 3430699
delta lambda=0.05
0.04-0.86

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 12.5000 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.2892 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 1.1840 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 50.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 2500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|-0.01168|4.35233|0.249377|2.26484|
|模拟|-0.00413|4.31998|0.186517|2.2212|
|误差|64.62%|0.74%|25.21%|1.93%|
### t235 13 3430821
delta lambda=0.05
0.04-0.86

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 8.7500 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.2420 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 1.6914 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 35.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 1750 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
**随着E_0减小，高频支由阻尼变为增长**

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|0.01402|3.88916|0.217628|2.08729|
|模拟|0.064356|3.80753|0.176653|2.06128|
|误差|359.04%|2.10%|18.83%|1.25%|
### t236 13 3430996
delta lambda=0.05
0.04-0.86

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 2.5000 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.1293 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 5.9200 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 10.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
**四支频率的振荡，如何理解？**
### t237 13 3431180
delta lambda=0.05
0.04-0.86

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 6.2500 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.2045 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 2.3680 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 25.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 1250 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |


|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.017505 | 3.56218 | 0.167516 | 1.89267 |
| 模拟 | 0.082605 | 3.52798 | 0.125764 | 1.88686 |
| 误差 | 371.90% | 0.96% | 24.92% | 0.31% |

### t238 13 
#### 33 3431341
delta lambda=0.05
0.04-0.86

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 8.2500 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.2350 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 1.7939 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 33.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 1650 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |


|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.014539 | 3.82448 | 0.209959 | 2.05527 |
| 模拟 | 0.057195 | 3.81403 | 0.150672 | 2.02687 |
| 误差 | 293.40% | 0.27% | 28.24% | 1.38% |

#### 31 3431345
delta lambda=0.05
0.04-0.86

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 7.7500 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.2277 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 1.9097 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 31.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 1550 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
#### 29 3431348
delta lambda=0.05
0.04-0.86

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 7.2500 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.2203 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 2.0414 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 29.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 1450 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
#### 27 3431350

^7ba67f

delta lambda=0.05
0.04-0.86

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 6.7500 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.2125 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 2.1926 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 27.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 1350 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
### t239 13 
delta lambda=0.05
0.04-0.86

#### 23 3432000
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 5.7500 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.1961 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 2.5739 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 23.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 1150 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
#### 21 3431806
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 5.2500 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.1874 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 2.8190 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 21.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 1050 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
#### 19 3431808
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 4.7500 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.1783 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 3.1158 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 19.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 950 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
#### 17 3431811
| 输入参数                     |         |          |         |  |           |          |          |           |
|--------------------------|---------|----------|---------|--|-----------|----------|----------|-----------|
| m_h                      | 1.00    | Lambda 0 | 0.1160  |  | 参量（理论对比）  |          |          |           |
| m_i                      | 1.00    |          |         |  | E_0/T_e   | 4.2500   | lambda_i | 0.0337    |
|                          |         | n_h/n_e  | 0.06    |  | E_c/T_e   | 14.8000  | lambda_h | 0.1686    |
| T_e/KeV                  | 4.00    |          |         |  |           |          |          |           |
| T_i/KeV                  | 0.60    | q        | 2.50    |  | tau_c     | 3.4824   |          |           |
|                          |         |          |         |  |           |          | nrx      | 191.1320  |
| E_0/KeV                  | 17.00   | e_i      | 1.00    |  | E_c/KeV   | 59.2000  |          |           |
|                          |         | e_h      | 1.00    |  |           |          |          |           |
| a                        | 0.60    |          |         |  | n_i/n_e   | 0.9400   |          |           |
| R                        | 3.75    |          |         |  |           |          |          |           |
| B0/T                     | 1.3750  |          |         |  |           |          |          |           |
|                          |         |          |         |  |           |          |          |           |
| perturbation wave length | 0.80    |          |         |  | 参量（程序调参）  |          |          |           |
|                          |         |          |         |  | T_h/KeV   | 850      |          |           |
|                          |         |          |         |  | v_0/C_s^h | 0.2      |
#### 15 3431814
| 输入参数                     |         |          |         |  |           |          |          |           |
|--------------------------|---------|----------|---------|--|-----------|----------|----------|-----------|
| m_h                      | 1.00    | Lambda 0 | 0.1160  |  | 参量（理论对比）  |          |          |           |
| m_i                      | 1.00    |          |         |  | E_0/T_e   | 3.7500   | lambda_i | 0.0337    |
|                          |         | n_h/n_e  | 0.06    |  | E_c/T_e   | 14.8000  | lambda_h | 0.1584    |
| T_e/KeV                  | 4.00    |          |         |  |           |          |          |           |
| T_i/KeV                  | 0.60    | q        | 2.50    |  | tau_c     | 3.9467   |          |           |
|                          |         |          |         |  |           |          | nrx      | 191.1320  |
| E_0/KeV                  | 15.00   | e_i      | 1.00    |  | E_c/KeV   | 59.2000  |          |           |
|                          |         | e_h      | 1.00    |  |           |          |          |           |
| a                        | 0.60    |          |         |  | n_i/n_e   | 0.9400   |          |           |
| R                        | 3.75    |          |         |  |           |          |          |           |
| B0/T                     | 1.3750  |          |         |  |           |          |          |           |
|                          |         |          |         |  |           |          |          |           |
| perturbation wave length | 0.80    |          |         |  | 参量（程序调参）  |          |          |           |
|                          |         |          |         |  | T_h/KeV   | 750      |          |           |
|                          |         |          |         |  | v_0/C_s^h | 0.2      |
#### 13 3431815
| 输入参数                     |         |          |         |  |           |          |          |           |
|--------------------------|---------|----------|---------|--|-----------|----------|----------|-----------|
| m_h                      | 1.00    | Lambda 0 | 0.1160  |  | 参量（理论对比）  |          |          |           |
| m_i                      | 1.00    |          |         |  | E_0/T_e   | 3.2500   | lambda_i | 0.0337    |
|                          |         | n_h/n_e  | 0.06    |  | E_c/T_e   | 14.8000  | lambda_h | 0.1475    |
| T_e/KeV                  | 4.00    |          |         |  |           |          |          |           |
| T_i/KeV                  | 0.60    | q        | 2.50    |  | tau_c     | 4.5538   |          |           |
|                          |         |          |         |  |           |          | nrx      | 191.1320  |
| E_0/KeV                  | 13.00   | e_i      | 1.00    |  | E_c/KeV   | 59.2000  |          |           |
|                          |         | e_h      | 1.00    |  |           |          |          |           |
| a                        | 0.60    |          |         |  | n_i/n_e   | 0.9400   |          |           |
| R                        | 3.75    |          |         |  |           |          |          |           |
| B0/T                     | 1.3750  |          |         |  |           |          |          |           |
|                          |         |          |         |  |           |          |          |           |
| perturbation wave length | 0.80    |          |         |  | 参量（程序调参）  |          |          |           |
|                          |         |          |         |  | T_h/KeV   | 650      |          |           |
|                          |         |          |         |  | v_0/C_s^h | 0.2      |

#### 11 3431817
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 2.7500 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.06 |  | E_c/T_e | 14.8000 | lambda_h | 0.1356 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 5.3818 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 11.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.9400 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 550 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
|  |  |  |  |  |  |  |  |  |
#### 修改浓度？？

### t240 13 3432362

^fededf

delta lambda=0.05
0.04-0.86
chen zhe文章中参数

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.3000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0367 |
|  |  | n_h/n_e | 0.20 |  | E_c/T_e | 14.8000 | lambda_h | 0.2000 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.3482 |  |  |
|  |  |  |  |  |  |  | nrx | 175.6145 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.8000 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 3.2620 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |


|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008755 | 3.66901 | 0.150686 | 1.25734 |
| 模拟 | 0.032723 | 3.60508 | 0.133024 | 1.25639 |
| 误差 | 273.75% | 1.74% | 11.72% | 0.08% |

### t241 13 [[#^7ba67f|E_0=27KeV(LHD)]] 3432363

delta lambda=0.05
0.04-0.86        n_h

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.1160 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 6.7500 | lambda_i | 0.0337 |
|  |  | n_h/n_e | 0.20 |  | E_c/T_e | 14.8000 | lambda_h | 0.2125 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 0.60 | q | 2.50 |  | tau_c | 2.1926 |  |  |
|  |  |  |  |  |  |  | nrx | 191.1320 |
| E_0/KeV | 27.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.8000 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 1350 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|0.054111|4.18434|0.301548|1.50807|
|模拟|0.126683|4.00067|0.307852|1.47053|
|误差|134.12%|4.39%|2.09%|2.49%|

### t242 13 [[#^fededf|t240]] 3432366

^2b8bd4

delta lambda=0.05
0.04-0.96 

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.3000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0368 |
|  |  | n_h/n_e | 0.20 |  | E_c/T_e | 14.8000 | lambda_h | 0.2009 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.3482 |  |  |
|  |  |  |  |  |  |  | nrx | 174.4104 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 1.26 |  |  |  | n_i/n_e | 0.8000 |  |  |
| R | 6.29 |  |  |  |  |  |  |  |
| B0/T | 1.3750 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.90 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
|  | γ1 | ω1 | γ2 | ω2 |  |  |  |  |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |  |  |  |  |
| 模拟 | 0.041755 | 3.58889 | 0.133384 | 1.28926 |  |  |  |  |
| 误差 | 373.08% | 2.17% | 11.48% | 2.54% |  |  |  |  |

### t243 E_0>>E_k 3432369
分布函数分布对比[[#^fededf|t240]]
delta lambda=0.05
0.04-0.86 

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.3000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 125.0000 | lambda_i | 0.0367 |
|  |  | n_h/n_e | 0.20 |  | E_c/T_e | 14.8000 | lambda_h | 0.3430 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.1184 |  |  |
|  |  |  |  |  |  |  | nrx | 175.6145 |
| E_0/KeV | 500.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.60 |  |  |  | n_i/n_e | 0.8000 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 3.2620 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 25000 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |

### t244 13 
delta lambda=0.05
0.04-0.86 
#### 170 3432973

^540dcb

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.3000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0366 |
|  |  | n_h/n_e | 0.20 |  | E_c/T_e | 14.8000 | lambda_h | 0.1997 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.3482 |  |  |
|  |  |  |  |  |  |  | nrx | 175.8657 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.98 |  |  |  | n_i/n_e | 0.8000 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 2.0000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
**chen zhe 文章较完美参数**


|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008732 | 3.66918 | 0.150689 | 1.25732 |
| 模拟 | 0.028561 | 3.61193 | 0.130148 | 1.25708 |
| 误差 | 227.09% | 1.56% | 13.63% | 0.02% |

#### 100 3432974
| 输入参数                     |         |          |         |  |           |          |          |           |
|--------------------------|---------|----------|---------|--|-----------|----------|----------|-----------|
| m_h                      | 1.00    | Lambda 0 | 0.3000  |  | 参量（理论对比）  |          |          |           |
| m_i                      | 1.00    |          |         |  | E_0/T_e   | 25.0000  | lambda_i | 0.0366    |
|                          |         | n_h/n_e  | 0.20    |  | E_c/T_e   | 14.8000  | lambda_h | 0.1532    |
| T_e/KeV                  | 4.00    |          |         |  |           |          |          |           |
| T_i/KeV                  | 4.00    | q        | 2.00    |  | tau_c     | 0.5920   |          |           |
|                          |         |          |         |  |           |          | nrx      | 175.8657  |
| E_0/KeV                  | 100.00  | e_i      | 1.00    |  | E_c/KeV   | 59.2000  |          |           |
|                          |         | e_h      | 1.00    |  |           |          |          |           |
| a                        | 0.98    |          |         |  | n_i/n_e   | 0.8000   |          |           |
| R                        | 3.75    |          |         |  |           |          |          |           |
| B0/T                     | 2.0000  |          |         |  |           |          |          |           |
|                          |         |          |         |  |           |          |          |           |
| perturbation wave length | 0.80    |          |         |  | 参量（程序调参）  |          |          |           |
|                          |         |          |         |  | T_h/KeV   | 5000     |          |           |
|                          |         |          |         |  | v_0/C_s^h | 0.2      |
#### 75 3433199
| 输入参数                     |         |          |         |  |           |          |          |           |
|--------------------------|---------|----------|---------|--|-----------|----------|----------|-----------|
| m_h                      | 1.00    | Lambda 0 | 0.3000  |  | 参量（理论对比）  |          |          |           |
| m_i                      | 1.00    |          |         |  | E_0/T_e   | 18.7500  | lambda_i | 0.0366    |
|                          |         | n_h/n_e  | 0.20    |  | E_c/T_e   | 14.8000  | lambda_h | 0.1327    |
| T_e/KeV                  | 4.00    |          |         |  |           |          |          |           |
| T_i/KeV                  | 4.00    | q        | 2.00    |  | tau_c     | 0.7893   |          |           |
|                          |         |          |         |  |           |          | nrx      | 175.8657  |
| E_0/KeV                  | 75.00   | e_i      | 1.00    |  | E_c/KeV   | 59.2000  |          |           |
|                          |         | e_h      | 1.00    |  |           |          |          |           |
| a                        | 0.98    |          |         |  | n_i/n_e   | 0.8000   |          |           |
| R                        | 3.75    |          |         |  |           |          |          |           |
| B0/T                     | 2.0000  |          |         |  |           |          |          |           |
|                          |         |          |         |  |           |          |          |           |
| perturbation wave length | 0.80    |          |         |  | 参量（程序调参）  |          |          |           |
|                          |         |          |         |  | T_h/KeV   | 3750     |          |           |
|                          |         |          |         |  | v_0/C_s^h | 0.2      |
#### 50 3433271
| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.3000 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 12.5000 | lambda_i | 0.0366   |
|                          |        | n_h/n_e  | 0.20   |     | E_c/T_e   | 14.8000 | lambda_h | 0.1083   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 1.1840  |          |          |
|                          |        |          |        |     |           |         | nrx      | 175.8657 |
| E_0/KeV                  | 50.00  | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 0.98   |          |        |     | n_i/n_e   | 0.8000  |          |          |
| R                        | 3.75   |          |        |     |           |         |          |          |
| B0/T                     | 2.0000 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.80   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 2500    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |


|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.007316 | 2.36454 | 0.078949 | 0.997318 |
| 模拟 | 0.032033 | 2.3396 | 0.082308 | 0.991386 |
| 误差 | 337.83% | 1.05% | 4.26% | 0.59% |

理论3:

| gamma3      | omega3  |
|-------------|---------|
| 3.61082E-05 | 3.01586 |


E_0=170kev为什么和理论结果不符合？什么效应影响？浓度？反环径比？
**可能是高阶的FOW效应使得170KeV仅有一支增长的模式   -->否定**
#### 170b 3435205
| 输入参数                     |         |          |         |  |           |          |          |           |
|--------------------------|---------|----------|---------|--|-----------|----------|----------|-----------|
| m_h                      | 1.00    | Lambda 0 | 0.3000  |  | 参量（理论对比）  |          |          |           |
| m_i                      | 1.00    |          |         |  | E_0/T_e   | 42.5000  | lambda_i | 0.0366    |
|                          |         | n_h/n_e  | 0.20    |  | E_c/T_e   | 14.8000  | lambda_h | 0.1997    |
| T_e/KeV                  | 4.00    |          |         |  |           |          |          |           |
| T_i/KeV                  | 4.00    | q        | 2.00    |  | tau_c     | 0.3482   |          |           |
|                          |         |          |         |  |           |          | nrx      | 175.8657  |
| E_0/KeV                  | 170.00  | e_i      | 1.00    |  | E_c/KeV   | 59.2000  |          |           |
|                          |         | e_h      | 1.00    |  |           |          |          |           |
| a                        | 0.98    |          |         |  | n_i/n_e   | 0.8000   |          |           |
| R                        | 6.75    |          |         |  |           |          |          |           |
| B0/T                     | 2.0000  |          |         |  |           |          |          |           |
|                          |         |          |         |  |           |          |          |           |
| perturbation wave length | 0.80    |          |         |  | 参量（程序调参）  |          |          |           |
|                          |         |          |         |  | T_h/KeV   | 8500     |          |           |
|                          |         |          |         |  | v_0/C_s^h | 0.2      |
1.25732 + 0.150689 I
3.66918 + 0.00873177 I
1.41242 - 0.763165 I
1.39391 - 0.670005 I

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|0.008732|3.66918|0.150689|1.25732|
|模拟|0.032832|3.60373|0.1328|1.25492|
|误差|276.00%|1.78%|11.87%|0.19%|

**反环径比减小对高频分支增长率误差影响很小，甚至为负效应**
#### 170c 3435242
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.3000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0366 |
|  |  | n_h/n_e | 0.40 |  | E_c/T_e | 14.8000 | lambda_h | 0.1997 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.3482 |  |  |
|  |  |  |  |  |  |  | nrx | 175.8657 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.98 |  |  |  | n_i/n_e | 0.6000 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 2.0000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|0.027031|4.20376|0.131429|1.01404|
|模拟|0.054694|4.03964|0.119789|0.988057|
|误差|102.34%|3.90%|8.86%|2.56%|
**高能离子浓度增大， 高频增长率、频率增大，低频分支均减小**
#### 170d  3435726
| 输入参数                     |         |          |         |  |           |           |          |           |
|--------------------------|---------|----------|---------|--|-----------|-----------|----------|-----------|
| m_h                      | 2.00    | Lambda 0 | 0.3000  |  | 参量（理论对比）  |           |          |           |
| m_i                      | 1.00    |          |         |  | E_0/T_e   | 42.5000   | lambda_i | 0.0260    |
|                          |         | n_h/n_e  | 0.40    |  | E_c/T_e   | 25.5085   | lambda_h | 0.2006    |
| T_e/KeV                  | 4.00    |          |         |  |           |           |          |           |
| T_i/KeV                  | 4.00    | q        | 2.00    |  | tau_c     | 0.6002    |          |           |
|                          |         |          |         |  |           |           | nrx      | 247.6477  |
| E_0/KeV                  | 170.00  | e_i      | 1.00    |  | E_c/KeV   | 102.0340  |          |           |
|                          |         | e_h      | 1.00    |  |           |           |          |           |
| a                        | 1.38    |          |         |  | n_i/n_e   | 0.6000    |          |           |
| R                        | 3.75    |          |         |  |           |           |          |           |
| B0/T                     | 2.0000  |          |         |  |           |           |          |           |
|                          |         |          |         |  |           |           |          |           |
| perturbation wave length | 0.80    |          |         |  | 参量（程序调参）  |           |          |           |
|                          |         |          |         |  | T_h/KeV   | 8500      |          |           |
|                          |         |          |         |  | v_0/C_s^h | 0.2       |
**模拟** 高能离子用氘非氕，高频分支转为阻尼
#### 170e 3436395
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 4.00 | Lambda 0 | 0.3000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0366 |
|  |  | n_h/n_e | 0.20 |  | E_c/T_e | 51.0170 | lambda_h | 0.1997 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 1.2004 |  |  |
|  |  |  |  |  |  |  | nrx | 175.8657 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 204.0681 |  |  |
|  |  | e_h | 2.00 |  |  |  |  |  |
| a | 0.98 |  |  |  | n_i/n_e | 0.6000 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 2.0000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
**理论**
1.37551 - 0.68815 I
1.30576 - 0.255321 
1.22984 + 0.227018 I
4.02794 - 0.030974 I

**模拟**
0.994
1.325

|  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- |
|  | γ1 | ω1 | γ2 | ω2 |  |
| 理论 |  |  |  |  |  |
| 模拟 | 0.034748 | 3.7024 | 0.100492 | 1.17719 |  |
| 误差 |  |  |  |  |  |

#### 170f 3436423
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.7000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0561 |
|  |  | n_h/n_e | 0.25 |  | E_c/T_e | 14.8000 | lambda_h | 0.2002 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.3482 |  |  |
|  |  |  |  |  |  |  | nrx | 114.8511 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.64 |  |  |  | n_i/n_e | 0.7500 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 2.0000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
**模拟**
频率 0.662V_ti/R0
#### 170g 3436432
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.5000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0436 |
|  |  | n_h/n_e | 0.25 |  | E_c/T_e | 14.8000 | lambda_h | 0.2011 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.3482 |  |  |
|  |  |  |  |  |  |  | nrx | 148.1159 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.94 |  |  |  | n_i/n_e | 0.7500 |  |  |
| R | 3.75 |  |  |  |  |  |  |  |
| B0/T | 2.0000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.70 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
**模拟**
频率 0.994V_ti/R0



170keV减小FOW效应：
#### 170h 3436444
| 输入参数                     |         |          |         |  |           |          |          |           |
|--------------------------|---------|----------|---------|--|-----------|----------|----------|-----------|
| m_h                      | 1.00    | Lambda 0 | 0.3000  |  | 参量（理论对比）  |          |          |           |
| m_i                      | 1.00    |          |         |  | E_0/T_e   | 42.5000  | lambda_i | 0.0198    |
|                          |         | n_h/n_e  | 0.25    |  | E_c/T_e   | 14.8000  | lambda_h | 0.1080    |
| T_e/KeV                  | 4.00    |          |         |  |           |          |          |           |
| T_i/KeV                  | 4.00    | q        | 2.00    |  | tau_c     | 0.3482   |          |           |
|                          |         |          |         |  |           |          | nrx      | 330.0218  |
| E_0/KeV                  | 170.00  | e_i      | 1.00    |  | E_c/KeV   | 59.2000  |          |           |
|                          |         | e_h      | 1.00    |  |           |          |          |           |
| a                        | 2.90    |          |         |  | n_i/n_e   | 0.7500   |          |           |
| R                        | 9.75    |          |         |  |           |          |          |           |
| B0/T                     | 2.0000  |          |         |  |           |          |          |           |
|                          |         |          |         |  |           |          |          |           |
| perturbation wave length | 0.50    |          |         |  | 参量（程序调参）  |          |          |           |
|                          |         |          |         |  | T_h/KeV   | 8500     |          |           |
|                          |         |          |         |  | v_0/C_s^h | 0.2      |
**理论**
1.18331 + 0.151521 I
1.26509 - 0.207138 I
3.87163 + 0.00408448 I

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|0.004085|3.87163|0.151521|1.18331|
|模拟|-0.02059|3.87669|0.125709|1.15228|
|误差|604.03%|0.13%|17.04%|2.62%|


扰动波长宽度增加到0.8，或减小到0.3
#### 170i 1513132

^3f82aa

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.3000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0184 |
|  |  | n_h/n_e | 0.25 |  | E_c/T_e | 14.8000 | lambda_h | 0.1004 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.3482 |  |  |
|  |  |  |  |  |  |  | nrx | 349.9369 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 1.95 |  |  |  | n_i/n_e | 0.7500 |  |  |
| R | 8.75 |  |  |  |  |  |  |  |
| B0/T | 2.0000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|-0.00355|3.87449|0.151512|1.18314|
|模拟|-0.02509|3.82546|0.12654|1.17773|
|误差|607.17%|1.27%|16.48%|0.46%|

#### 170j 1513262
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.3000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0192 |
|  |  | n_h/n_e | 0.25 |  | E_c/T_e | 14.8000 | lambda_h | 0.1046 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.3482 |  |  |
|  |  |  |  |  |  |  | nrx | 349.4555 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 4.99 |  |  |  | n_i/n_e | 0.7500 |  |  |
| R | 24.00 |  |  |  |  |  |  |  |
| B0/T | 2.0000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.30 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
R_R = 0.35
L_R = 0.05

|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|-0.00384|3.87294|0.151517|1.18323|
|模拟|-0.01354|3.8477|0.135682|1.15664|
|误差|252.56%|0.65%|10.45%|2.25%|


#### 170k 1513311
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.3000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0192 |
|  |  | n_h/n_e | 0.25 |  | E_c/T_e | 14.8000 | lambda_h | 0.1046 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.3482 |  |  |
|  |  |  |  |  |  |  | nrx | 349.4555 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 4.99 |  |  |  | n_i/n_e | 0.7500 |  |  |
| R | 24.00 |  |  |  |  |  |  |  |
| B0/T | 2.0000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.30 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
R_R = 0.65
 L_R = 0.35
 **理论**
 1.18323 + 0.151517 I
 3.87294 + 0.00384042 I
##### 170l 1525489
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.3000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0184 |
|  |  | n_h/n_e | 0.25 |  | E_c/T_e | 14.8000 | lambda_h | 0.1004 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.3482 |  |  |
|  |  |  |  |  |  |  | nrx | 363.9869 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.99 |  |  |  | n_i/n_e | 0.7500 |  |  |
| R | 4.90 |  |  |  |  |  |  |  |
| B0/T | 10.5000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.30 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |

**理论**
3.87449 - 0.00354847 I
3.87449 + 0.0035482 I
1.18314 + 0.151512 I


|   |   |   |   |   |
|---|---|---|---|---|
||γ1|ω1|γ2|ω2|
|理论|-0.00355|3.87449|0.151512|1.18314|
|模拟|-0.01221|3.8199|0.104441|1.19702|
|误差|244.09%|1.41%|31.07%|1.17%|


**初始能量越高，高频模的阻尼效应越强**
##### 170m 1525548
| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.3000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0197 |
|  |  | n_h/n_e | 0.25 |  | E_c/T_e | 14.8000 | lambda_h | 0.1076 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.3482 |  |  |
|  |  |  |  |  |  |  | nrx | 324.9883 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 1.50 |  |  |  | n_i/n_e | 0.7500 |  |  |
| R | 7.20 |  |  |  |  |  |  |  |
| B0/T | 2.0000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.97 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
R_R = 0.985
 L_R = 0.015

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | -0.00406 | 3.87179 | 0.15152 | 1.1833 |
| 模拟 | -0.02277 | 3.80359 | 0.124209 | 1.18724 |
| 误差 | 461.32% | 1.76% | 18.02% | 0.33% |

#### 170n 1525635 [[#^3f82aa|170i]]
R_R =0.90
 L_R =0.10

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.3000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 42.5000 | lambda_i | 0.0184 |
|  |  | n_h/n_e | 0.25 |  | E_c/T_e | 14.8000 | lambda_h | 0.1004 |
| T_e/KeV | 4.00 |  |  |  |  |  |  |  |
| T_i/KeV | 4.00 | q | 2.00 |  | tau_c | 0.3482 |  |  |
|  |  |  |  |  |  |  | nrx | 349.9369 |
| E_0/KeV | 170.00 | e_i | 1.00 |  | E_c/KeV | 59.2000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 1.95 |  |  |  | n_i/n_e | 0.7500 |  |  |
| R | 8.75 |  |  |  |  |  |  |  |
| B0/T | 2.0000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.80 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 8500 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |

### t245 13 [[#^2b8bd4|t242]] 1533656

^0f8be0

| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.3000 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 42.5000 | lambda_i | 0.0368   |
|                          |        | n_h/n_e  | 0.20   |     | E_c/T_e   | 14.8000 | lambda_h | 0.2009   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.3482  |          |          |
|                          |        |          |        |     |           |         | nrx      | 174.4104 |
| E_0/KeV                  | 170.00 | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 1.26   |          |        |     | n_i/n_e   | 0.8000  |          |          |
| R                        | 6.29   |          |        |     |           |         |          |          |
| B0/T                     | 1.3750 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.90   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 8500    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |
nrx = 180
ntm = 1
nalpha = 2
ntheta = 16
nvpara = 160
muType = 3
nmu = 32
0.04-0.96 
delta lambda=0.04


|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.046078 | 3.5959 | 0.136511 | 1.26636 |
| 误差 | 422.07% | 1.98% | 9.40% | 0.71% |

#### delta lambda=0.02 1534037

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.090935 | 3.58363 | 0.153323 | 1.26868 |
| 误差 | 930.31% | 2.31% | 1.75% | 0.90% |
#### delta lambda=0.03 1534042

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.063937 | 3.59253 | 0.144458 | 1.26779 |
| 误差 | 624.42% | 2.07% | 4.13% | 0.83% |
#### delta lambda=0.055 1538807


|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.034067 | 3.59926 | 0.131029 | 1.25842 |
| 误差 | 285.98% | 1.89% | 13.04% | 0.08% |

#### delta lambda=0.06  1538815


|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.03234 | 3.60042 | 0.13107 | 1.25542 |
| 误差 | 266.41% | 1.86% | 13.01% | 0.16% |

#### delta lambda=0.1 1540094

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.028731 | 3.60397 | 0.135947 | 1.2448 |
| 误差 | 225.53% | 1.76% | 9.78% | 1.00% |
#### delta lambda=0.15 1540918

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.025752 | 3.62346 | 0.14509 | 1.23244 |
| 误差 | 191.78% | 1.23% | 3.71% | 1.98% |

#### delta lambda=0.2 1540920

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.015112 | 3.62171 | 0.149053 | 1.20936 |
| 误差 | 71.22% | 1.28% | 1.08% | 3.82% |

**是否因为理论值对应的delta lambda偏移**

#### delta lambda=0.3 1541584

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | -0.03694 | 3.61022 | 0.109847 | 1.12594 |
| 误差 | 518.55% | 1.59% | 27.10% | 10.45% |

#### delta lambda=0.4 1541591

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | -0.02701 | 3.51607 | 0.029713 | 1.07591 |
| 误差 | 406.06% | 4.16% | 80.28% | 14.43% |
#### simTotalTime = 20.1 maxStep = 50000

#### delta lambda=0.21 1541898

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.012089 | 3.61995 | 0.148827 | 1.20309 |
| 误差 | 36.97% | 1.32% | 1.23% | 4.32% |
#### delta lambda=0.25 1541895

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | -0.00142 | 3.60922 | 0.141219 | 1.17397 |
| 误差 | 116.10% | 1.62% | 6.28% | 6.63% |
#### delta lambda=0.26 1541896
|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | -0.00491 | 3.60573 | 0.137358 | 1.16607 |
| 误差 | 155.64% | 1.71% | 8.84% | 7.26% |
#### delta lambda=0.01 1541925
|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | -0.02994 | 3.77688 | 0.110956 | 1.34072 |
| 误差 | 439.21% | 2.95% | 26.36% | 6.63% |
#### delta lambda=0.22 1542198

|     | γ1       | ω1      | γ2       | ω2      |
| --- | -------- | ------- | -------- | ------- |
| 理论  | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟  | 0.008873 | 3.61779 | 0.148007 | 1.19634 |
| 误差  | 0.53%    | 1.38%   | 1.77%    | 4.85%   |

### t246 13 [[#^540dcb|t244 chen zhe]] 1540935
| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.3000 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 42.5000 | lambda_i | 0.0366   |
|                          |        | n_h/n_e  | 0.20   |     | E_c/T_e   | 14.8000 | lambda_h | 0.1997   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.3482  |          |          |
|                          |        |          |        |     |           |         | nrx      | 175.8657 |
| E_0/KeV                  | 170.00 | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 0.98   |          |        |     | n_i/n_e   | 0.8000  |          |          |
| R                        | 3.75   |          |        |     |           |         |          |          |
| B0/T                     | 2.0000 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.80   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 8500    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |
delta lambda=0.04
0.04-0.86

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008732 | 3.66918 | 0.150689 | 1.25732 |
| 模拟 | 0.037722 | 3.60873 | 0.134979 | 1.26275 |
| 误差 | 331.99% | 1.65% | 10.43% | 0.43% |
#### delta lambda=0.03 1541899

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008732 | 3.66918 | 0.150689 | 1.25732 |
| 模拟 | 0.052894 | 3.60359 | 0.143359 | 1.26494 |
| 误差 | 505.75% | 1.79% | 4.86% | 0.61% |
#### delta lambda=0.02 1541901
|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008732 | 3.66918 | 0.150689 | 1.25732 |
| 模拟 | 0.037722 | 3.60873 | 0.134979 | 1.26275 |
| 误差 | 331.99% | 1.65% | 10.43% | 0.43% |
#### delta lambda=0.01 1541922
|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008732 | 3.66918 | 0.150689 | 1.25732 |
| 模拟 | -0.03093 | 3.78936 | 0.113452 | 1.31306 |
| 误差 | 454.26% | 3.28% | 24.71% | 4.43% |

### t247 13 1543198
[[chen2020verification]]
delta lambda=0.05
lambda_0=0.5
0.19-0.81
ntheta=16 nvpara 160  nmu  32

| 输入参数 |  |  |  |  |  |  |  |  |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| m_h | 1.00 | Lambda 0 | 0.5000 |  | 参量（理论对比） |  |  |  |
| m_i | 1.00 |  |  |  | E_0/T_e | 63.3600 | lambda_i | 0.0469 |
|  |  | n_h/n_e | 0.08 |  | E_c/T_e | 14.8000 | lambda_h | 0.2640 |
| T_e/KeV | 1.00 |  |  |  |  |  |  |  |
| T_i/KeV | 1.00 | q | 4.00 |  | tau_c | 0.2336 |  |  |
|  |  |  |  |  |  |  | nrx | 138.3991 |
| E_0/KeV | 63.36 | e_i | 1.00 |  | E_c/KeV | 14.8000 |  |  |
|  |  | e_h | 1.00 |  |  |  |  |  |
| a | 0.51 |  |  |  | n_i/n_e | 0.9200 |  |  |
| R | 1.70 |  |  |  |  |  |  |  |
| B0/T | 2.0000 |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |
| perturbation wave length | 0.60 |  |  |  | 参量（程序调参） |  |  |  |
|  |  |  |  |  | T_h/KeV | 3168 |  |  |
|  |  |  |  |  | v_0/C_s^h | 0.2 |  |  |
**theory**
{\[Omega] -> 0.930694 + 0.0797932 I}
{\[Omega] -> 1.93906 + 0.142015 I}
#### ntheta=64 1543426
结果变化明显
#### GAM 1543434
**theory** 1.70162 -0.00030778 I
0.5 a: omega=1.67092，gamma=-0.00517669
$\omega_G=\sqrt{(7\times T_i/2+2\times T_e)/(m_i\times R_0^2)}=1.6583$

## 2024-02-22
#### ntheta=1320 1550037
not run
#### nmu=128 ntheta=64 1554609
继续增加速度网格



### t248 13 [[#^0f8be0|t245]] 
速度网格调节，原参数：ntheta 16   nvpara 160  nmu  32
Δλ=0.05
#### 1.仅调节ntheta
##### ntheta=64 1549703
**theory**:frequency:1.325;3.312
#### 2.仅调节nmu
##### nmu=96 1549761
**有效减小误差** mu网格影响很大

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.020025 | 3.5855 | 0.140233 | 1.2522 |
| 误差 | 126.89% | 2.26% | 6.93% | 0.41% |
###### time max=40.2
nmu=32 [[#^2b8bd4|t242]]
nmu=96 1554720
##### nmu=128 1554614
**nmu=96**fit well

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.020362 | 3.58664 | 0.139627 | 1.25233 |
| 误差 | 130.71% | 2.23% | 7.33% | 0.40% |

#### 3.仅调节nvpara
##### nvpara=256 1549770
区别不大

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.038274 | 3.61298 | 0.128798 | 1.26149 |
| 误差 | 333.65% | 1.51% | 14.52% | 0.33% |
|  |  |  |  |  |
|  | γ1 | ω1 | γ2 | ω2 |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.044757 | 3.60058 | 0.130647 | 1.29107 |
| 误差 | 407.10% | 1.85% | 13.29% | 2.68% |
## 2024-02-23
#### 2nmu 
**nmu=96**网格对Δλ>0.05是否足够描述分布函数变化?
##### Δλ=0.1 1554677
网格加大，fit结果有极小区别，**认为模拟结果相同**

|  | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.028664 | 3.60306 | 0.137246 | 1.2462 |
| 误差 | 224.77% | 1.78% | 8.91% | 0.89% |
##### Δλ=0.2 1554706
| cr | γ1        | ω1      | γ2       | ω2      |
|----|-----------|---------|----------|---------|
| 理论 | 0.008826  | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.0158597 | 3.61416 | 0.149225 | 1.2092  |
| 误差 | 79.69%    | 1.48%   | 0.96%    | 3.83%   |

**nmu=96**网格对Δλ<0.05是否足够描述分布函数变化?
##### aΔλ=0.01 1554710
1-70 [0196](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct248%5Cnmulambda%5Cb%5C0196)

| cr  | γ1         | ω1      | γ2       | ω2      |
| --- | ---------- | ------- | -------- | ------- |
| 理论  | 0.008826   | 3.66852 | 0.150679 | 1.25738 |
| 模拟  | 0.00844335 | 3.63024 | 0.131256 | 1.24836 |
| 误差  | 4.34%      | 1.04%   | 12.89%   | 0.72%   |
|     |            |         |          |         |
1-100

| cr | γ1        | ω1      | γ2       | ω2      |
|----|-----------|---------|----------|---------|
| 理论 | 0.008826  | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.0300528 | 3.60211 | 0.158044 | 1.24354 |
| 误差 | 240.50%   | 1.81%   | 4.89%    | 1.10%   |

##### bΔλ=0.01 1554713
**nmu=128**
[01128](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct248%5Cnmulambda%5Cb%5C01128)
1-110

| cr | γ1 | ω1 | γ2 | ω2 |
| ---- | ---- | ---- | ---- | ---- |
| 理论 | 0.008826 | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.0295832 | 3.57469 | 0.136906 | 1.24338 |
| 误差 | 235.18% | 2.56% | 9.14% | 1.11% |


1-70

| cr  | γ1        | ω1      | γ2       | ω2      |
| --- | --------- | ------- | -------- | ------- |
| 理论  | 0.008826  | 3.66852 | 0.150679 | 1.25738 |
| 模拟  | 0.0106694 | 3.60582 | 0.125049 | 1.23667 |
| 误差  | 20.89%    | 1.71%   | 17.01%   | 1.65%   |


**拟合选择fit较好的最大范围** 利于比较
低频未增长起来的范围内拟合高频阻尼率拟合结果较准确，低频占主导时高频分支拟合不准确

验证模拟结果可靠：
### t249 13 
[[#bΔλ=0.01 1554713]]

| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.3000 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 42.5000 | lambda_i | 0.0368   |
|                          |        | n_h/n_e  | 0.20   |     | E_c/T_e   | 14.8000 | lambda_h | 0.2009   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.3482  |          |          |
|                          |        |          |        |     |           |         | nrx      | 174.4104 |
| E_0/KeV                  | 170.00 | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 1.26   |          |        |     | n_i/n_e   | 0.8000  |          |          |
| R                        | 6.29   |          |        |     |           |         |          |          |
| B0/T                     | 1.3750 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.90   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 8500    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |
邻近Δλ的模拟
#### Δλ=0.011 1593475
[deltalam011](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct249%5Cdeltalam011)
1-150

| cr | γ1        | ω1      | γ2       | ω2      |
|----|-----------|---------|----------|---------|
| 理论 | 0.008826  | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.0233046 | 3.56933 | 0.140326 | 1.25224 |
| 误差 | 164.04%   | 2.70%   | 6.87%    | 0.41%   |
1-110

| cr | γ1        | ω1      | γ2       | ω2      |
|----|-----------|---------|----------|---------|
| 理论 | 0.008826  | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.0288053 | 3.5805  | 0.135873 | 1.2495  |
| 误差 | 226.37%   | 2.40%   | 9.83%    | 0.63%   |
1-70

| cr | γ1         | ω1      | γ2       | ω2      |
|----|------------|---------|----------|---------|
| 理论 | 0.008826   | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.00930597 | 3.61076 | 0.122108 | 1.24404 |
| 误差 | 5.44%      | 1.57%   | 18.96%   | 1.06%   |

#### Δλ=0.012 1594158
[deltalam012](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct249%5Cdeltalam012)
1-110

| cr | γ1        | ω1      | γ2       | ω2      |
|----|-----------|---------|----------|---------|
| 理论 | 0.008826  | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.0285406 | 3.58422 | 0.135904 | 1.25339 |
| 误差 | 223.37%   | 2.30%   | 9.81%    | 0.32%   |
1-70

| cr | γ1         | ω1      | γ2       | ω2      |
|----|------------|---------|----------|---------|
| 理论 | 0.008826   | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.00858085 | 3.61415 | 0.120711 | 1.24894 |
| 误差 | 2.78%      | 1.48%   | 19.89%   | 0.67%   |

#### Δλ=0.013 1594314
[deltalam013](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct249%5Cdeltalam013)
1-110

| cr  | γ1        | ω1      | γ2       | ω2      |
| --- | --------- | ------- | -------- | ------- |
| 理论  | 0.008826  | 3.66852 | 0.150679 | 1.25738 |
| 模拟  | 0.0285133 | 3.58664 | 0.136151 | 1.25576 |
| 误差  | 223.06%   | 2.23%   | 9.64%    | 0.13%   |

1-70

| cr | γ1         | ω1      | γ2       | ω2      |
|----|------------|---------|----------|---------|
| 理论 | 0.008826   | 3.66852 | 0.150679 | 1.25738 |
| 模拟 | 0.00820019 | 3.61653 | 0.119976 | 1.252   |
| 误差 | 7.09%      | 1.42%   | 20.38%   | 0.43%   |

#### 改变E_0，lambda0 ，B
##### 1 1595134
| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.3000 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 30.0000 | lambda_i | 0.0368   |
|                          |        | n_h/n_e  | 0.20   |     | E_c/T_e   | 14.8000 | lambda_h | 0.1688   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.4933  |          |          |
|                          |        |          |        |     |           |         | nrx      | 174.4104 |
| E_0/KeV                  | 120.00 | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 1.26   |          |        |     | n_i/n_e   | 0.8000  |          |          |
| R                        | 6.29   |          |        |     |           |         |          |          |
| B0/T                     | 1.3750 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.90   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 6000    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |
[1](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct249%5C1)
1-110

| cr | γ1         | ω1      | γ2       | ω2      |
|----|------------|---------|----------|---------|
| 理论 | 0.00850781 | 3.20993 | 0.14411  | 1.19837 |
| 模拟 | 0.0175805  | 3.16295 | 0.127652 | 1.19609 |
| 误差 | 106.64%    | 1.46%   | 11.42%   | 0.19%   |

1-50

| cr | γ1         | ω1      | γ2      | ω2      |
|----|------------|---------|---------|---------|
| 理论 | 0.00850781 | 3.20993 | 0.14411 | 1.19837 |
| 模拟 | -0.0137486 | 3.16643 | 0.06043 | 1.19145 |
| 误差 | 261.60%    | 1.36%   | 58.07%  | 0.58%   |
1-100

| cr  | γ1         | ω1      | γ2       | ω2      |
| --- | ---------- | ------- | -------- | ------- |
| 理论  | 0.00850781 | 3.20993 | 0.14411  | 1.19837 |
| 模拟  | 0.0088873  | 3.16306 | 0.128995 | 1.2049  |
| 误差  | 4.46%      | 1.46%   | 10.49%   | 0.54%   |


##### 2 1595651
| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.5000 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 30.0000 | lambda_i | 0.0368   |
|                          |        | n_h/n_e  | 0.20   |     | E_c/T_e   | 14.8000 | lambda_h | 0.1426   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.4933  |          |          |
|                          |        |          |        |     |           |         | nrx      | 174.4104 |
| E_0/KeV                  | 120.00 | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 1.26   |          |        |     | n_i/n_e   | 0.8000  |          |          |
| R                        | 6.29   |          |        |     |           |         |          |          |
| B0/T                     | 1.3750 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.90   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 6000    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |
[2](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct249%5C2)
1-110

| cr | γ1        | ω1      | γ2       | ω2      |
|----|-----------|---------|----------|---------|
| 理论 | 0.0259859 | 2.8929  | 0.195744 | 1.04817 |
| 模拟 | 0.082781  | 2.86283 | 0.166533 | 1.04474 |
| 误差 | 218.56%   | 1.04%   | 14.92%   | 0.33%   |

1-70

| cr | γ1        | ω1      | γ2       | ω2      |
|----|-----------|---------|----------|---------|
| 理论 | 0.0259859 | 2.8929  | 0.195744 | 1.04817 |
| 模拟 | 0.0356669 | 2.85923 | 0.159072 | 1.05105 |
| 误差 | 37.25%    | 1.16%   | 18.73%   | 0.27%   |

1-150

| cr | γ1        | ω1      | γ2       | ω2      |
|----|-----------|---------|----------|---------|
| 理论 | 0.0259859 | 2.8929  | 0.195744 | 1.04817 |
| 模拟 | 0.0892875 | 2.81373 | 0.162387 | 1.03875 |
| 误差 | 243.60%   | 2.74%   | 17.04%   | 0.90%   |

1-160

| cr | γ1        | ω1     | γ2       | ω2      |
|----|-----------|--------|----------|---------|
| 理论 | 0.0259859 | 2.8929 | 0.195744 | 1.04817 |
| 模拟 | 0.0857464 | 2.805  | 0.169877 | 1.03823 |
| 误差 | 229.97%   | 3.04%  | 13.21%   | 0.95%   |

##### 3 1596291
[3](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct249%5C3)

| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.5000 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 30.0000 | lambda_i | 0.0506   |
|                          |        | n_h/n_e  | 0.20   |     | E_c/T_e   | 14.8000 | lambda_h | 0.1961   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.4933  |          |          |
|                          |        |          |        |     |           |         | nrx      | 126.8439 |
| E_0/KeV                  | 120.00 | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 1.26   |          |        |     | n_i/n_e   | 0.8000  |          |          |
| R                        | 6.29   |          |        |     |           |         |          |          |
| B0/T                     | 1.0000 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.90   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 6000    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |
1-110

| cr | γ1        | ω1      | γ2       | ω2      |
|----|-----------|---------|----------|---------|
| 理论 | 0.0443095 | 2.85781 | 0.193174 | 1.05179 |
| 模拟 | 0.114954  | 2.72669 | 0.163952 | 1.04145 |
| 误差 | 159.43%   | 4.59%   | 15.13%   | 0.98%   |

1-60

| cr | γ1        | ω1      | γ2       | ω2      |
|----|-----------|---------|----------|---------|
| 理论 | 0.0443095 | 2.85781 | 0.193174 | 1.05179 |
| 模拟 | 0.0762223 | 2.78343 | 0.117303 | 1.05589 |
| 误差 | 72.02%    | 2.60%   | 39.28%   | 0.39%   |

1-200

| cr | γ1        | ω1      | γ2       | ω2      |
|----|-----------|---------|----------|---------|
| 理论 | 0.0443095 | 2.85781 | 0.193174 | 1.05179 |
| 模拟 | 0.0492269 | 2.73846 | 0.167409 | 1.0466  |
| 误差 | 11.10%    | 4.18%   | 13.34%   | 0.49%   |

## 2024-02-26

### t250 13 [[#t245 13 2b8bd4 t242 1533656|t245]]
num=128
[t250](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct250)

| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.3000 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 42.5000 | lambda_i | 0.0368   |
|                          |        | n_h/n_e  | 0.20   |     | E_c/T_e   | 14.8000 | lambda_h | 0.2009   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.3482  |          |          |
|                          |        |          |        |     |           |         | nrx      | 174.4104 |
| E_0/KeV                  | 170.00 | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 1.26   |          |        |     | n_i/n_e   | 0.8000  |          |          |
| R                        | 6.29   |          |        |     |           |         |          |          |
| B0/T                     | 1.3750 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.90   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 8500    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |

### t251 13 [[#170 3432973|chenzhe]] 
nmu=128       0.04-0.86 
[t251](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct251)

| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.3000 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 42.5000 | lambda_i | 0.0366   |
|                          |        | n_h/n_e  | 0.20   |     | E_c/T_e   | 14.8000 | lambda_h | 0.1997   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.3482  |          |          |
|                          |        |          |        |     |           |         | nrx      | 175.8657 |
| E_0/KeV                  | 170.00 | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 0.98   |          |        |     | n_i/n_e   | 0.8000  |          |          |
| R                        | 3.75   |          |        |     |           |         |          |          |
| B0/T                     | 2.0000 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.80   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 8500    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |
增长率为负值认为是增长率太小拟合出错导致，改变浓度可增大高频分支增长率
## 2024-02-29
### t252 13 1672395
deltalambda=0.01
[t252](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct252)

| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.3000 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 42.5000 | lambda_i | 0.0366   |
|                          |        | n_h/n_e  | 0.50   |     | E_c/T_e   | 14.8000 | lambda_h | 0.1997   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.3482  |          |          |
|                          |        |          |        |     |           |         | nrx      | 175.8657 |
| E_0/KeV                  | 170.00 | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 0.98   |          |        |     | n_i/n_e   | 0.5000  |          |          |
| R                        | 3.75   |          |        |     |           |         |          |          |
| B0/T                     | 2.0000 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.80   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 8500    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |
### t253 13 1672826
deltalambda=0.01
[t253](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct253)

| 输入参数                     |        |          |        |     |           |          |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | -------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.3000 |     | 参量（理论对比）  |          |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 425.0000 | lambda_i | 0.0366   |
|                          |        | n_h/n_e  | 0.50   |     | E_c/T_e   | 14.8000  | lambda_h | 0.1997   |
| T_e/KeV                  | 0.40   |          |        |     |           |          |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.0348   |          |          |
|                          |        |          |        |     |           |          | nrx      | 175.8657 |
| E_0/KeV                  | 170.00 | e_i      | 1.00   |     | E_c/KeV   | 5.9200   |          |          |
|                          |        | e_h      | 1.00   |     |           |          |          |          |
| a                        | 0.98   |          |        |     | n_i/n_e   | 0.5000   |          |          |
| R                        | 3.75   |          |        |     |           |          |          |          |
| B0/T                     | 2.0000 |          |        |     |           |          |          |          |
|                          |        |          |        |     |           |          |          |          |
| perturbation wave length | 0.80   |          |        |     | 参量（程序调参）  |          |          |          |
|                          |        |          |        |     | T_h/KeV   | 8500     |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2      |          |          |
#### flat 0density0  1685130
[flat](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct253%5Cflat)
#### not flat thermal ion density  1685136
频率基本不变，高频分支增长率略有下降[nonflat](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct253%5Cnonflat)

增加本底等体温度梯度 1708078
无高能离子情况  1708127
无温度梯度情况 对比  1719294
&Par
  !Electron
  mess   = 5.446d-4  ! 1/1836.152
  charge = -1.0
  vMax   = 4.0
  PType  = 4
  Den    = 1.0
  Upa    = 0.0
  Tem    = 0.3654
/

&Par
  !Ion1
  mess   = 1.0
  charge = 1.0
  vMax   = 4.0
  PType  = 4
  Den    = 0.5
  Upa    = 0.0
  Tem    = 3.6543
/
          
&Par     
  !Ion2   
  mess   = 1.0
  charge = 1.0
  vMax   = 0.25
  PType  = 6 
  Den    = 0.5
  Upa    = 0.2
  Tem    = 7765.5

温度梯度 kappaT = 9.96_dp 1733658
[nonflattemp2](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct253%5Cnonflattemp2)

## 2024-03-07

q   c0=2  c1=1   1785167[nablaq](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct253%5Cnablaq)

n0 T0 q GAM  1799113 [nablaqGAM](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct253%5CnablaqGAM)

### t254 13 1767272
Δλ=0.2
调高高频阻尼后扫Δλ

| 输入参数                     |        |          |        |     |           |          |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | -------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.3000 |     | 参量（理论对比）  |          |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 425.0000 | lambda_i | 0.0366   |
|                          |        | n_h/n_e  | 0.50   |     | E_c/T_e   | 14.8000  | lambda_h | 0.1997   |
| T_e/KeV                  | 0.40   |          |        |     |           |          |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.0348   |          |          |
|                          |        |          |        |     |           |          | nrx      | 175.8657 |
| E_0/KeV                  | 170.00 | e_i      | 1.00   |     | E_c/KeV   | 5.9200   |          |          |
|                          |        | e_h      | 1.00   |     |           |          |          |          |
| a                        | 0.98   |          |        |     | n_i/n_e   | 0.5000   |          |          |
| R                        | 3.75   |          |        |     |           |          |          |          |
| B0/T                     | 2.0000 |          |        |     |           |          |          |          |
|                          |        |          |        |     |           |          |          |          |
| perturbation wave length | 0.80   |          |        |     | 参量（程序调参）  |          |          |          |
|                          |        |          |        |     | T_h/KeV   | 8500     |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2      |          |          |

Δλ=0.17  1773671

## 2024-03-13

### t255 13 1961928
0.02-0.92
Δλ=0.05

| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.3000 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 42.5000 | lambda_i | 0.0371   |
|                          |        | n_h/n_e  | 0.25   |     | E_c/T_e   | 14.8000 | lambda_h | 0.2025   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.3482  |          |          |
|                          |        |          |        |     |           |         | nrx      | 173.0262 |
| E_0/KeV                  | 170.00 | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 1.25   |          |        |     | n_i/n_e   | 0.7500  |          |          |
| R                        | 3.75   |          |        |     |           |         |          |          |
| B0/T                     | 1.3750 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.90   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 8500    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |
## 2024-03-14

Te=0.4 1971181

E0=75kev Te=0.4 1971427

### t256 13  [[#t245 13 2b8bd4 t242 1533656|t245]]
Δλ=0.05 1972178
Δλ=0.2 1972187

GAM 1978983
## 2024-03-18
nonlinear   delta λ=0.05  2007527
 nalpha=2-> nalpha = 382


### t257 13 2129747
[[#t251 13 170 3432973 chenzhe|t251]]
Δλ=0.05
ntheta=16
[t257](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct257)

| cr  | ω1      | γ1      | ω2      | γ2      |
| --- | ------- | ------- | ------- | ------- |
| 理论  | 1.25732 | 0.15069 | 3.66918 | 0.00873 |
| 模拟  | 1.23928 | 0.14696 | 3.60356 | 0.02724 |
| 误差  | 1.43%   | 2.47%   | 1.79%   | 211.95% |

### t258 13 2131074
[[#t251 13 170 3432973 chenzhe|t251]]
Δλ=0.05
ntheta=64[t258](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct258)

| cr | ω1       | γ1       | ω2       | γ2       |
|----|----------|----------|----------|----------|
| 理论 | 1.25732  | 0.15069  | 3.66918  | 0.00873  |
| 模拟 | 1.26224  | 0.15019  | 3.65739  | 0.02440  |
| 误差 | 0.39%    | 0.33%    | 0.32%    | 179.49%  |

### t259 13 2131113

| 输入参数                     |          |          |        |     |          |           |          |          |
| ------------------------ | -------- | -------- | ------ | --- | -------- | --------- | -------- | -------- |
| m_h                      | 1.00     | Lambda 0 | 0.3000 |     | 参量（理论对比） |           |          |          |
| m_i                      | 1.00     |          |        |     | E_0/T_e  | 1700.0000 | lambda_i | 0.0366   |
|                          |          | n_h/n_e  | 0.20   |     | E_c/T_e  | 14.8000   | lambda_h | 0.1997   |
| ==T_e/KeV==              | ==0.10== |          |        |     |          |           |          |          |
| T_i/KeV                  | 4.00     | q        | 2.00   |     | tau_c    | 0.0087    |          |          |
|                          |          |          |        |     |          |           | nrx      | 175.8657 |
| E_0/KeV                  | 170.00   | e_i      | 1.00   |     | E_c/KeV  | 1.4800    |          |          |
|                          |          | e_h      | 1.00   |     |          |           |          |          |
| a                        | 0.98     |          |        |     | n_i/n_e  | 0.8000    |          |          |
| R                        | 3.75     |          |        |     |          |           |          |          |
| B0/T                     | 2.0000   |          |        |     |          |           |          |          |
|                          |          |          |        |     |          |           |          |          |
| perturbation wave length | **0.80** |          |        |     | 参量（程序调参） |           |          |          |
|                          |          |          |        |     | T_h/KeV  | 8500      |          |          |
|                          |          |          |        |     |          |           |          |          |
$E_c\ll E_0$，$\Lambda0$较小出现阻尼；温度简化理论对比无意义
Δλ=0.1
ntheta=16

## 2024-03-24
### t260 13 2143658
[[#t251 13 170 3432973 chenzhe|t251]]
仅改变lambda=0.75观察残余流
delta lambda=0.05

| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.7500 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 42.5000 | lambda_i | 0.0366   |
|                          |        | n_h/n_e  | 0.20   |     | E_c/T_e   | 14.8000 | lambda_h | 0.1194   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.3482  |          |          |
|                          |        |          |        |     |           |         | nrx      | 175.8657 |
| E_0/KeV                  | 170.00 | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 0.98   |          |        |     | n_i/n_e   | 0.8000  |          |          |
| R                        | 3.75   |          |        |     |           |         |          |          |
| B0/T                     | 2.0000 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.80   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 8500    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |
**扰动在边缘被激发？时间增加测试**
### t261 13 NLT 2144228

| 输入参数                     |        |          |        |     |           |         |          |          |
| ------------------------ | ------ | -------- | ------ | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0 | 0.3000 |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |          |        |     | E_0/T_e   | 42.5000 | lambda_i | 0.0368   |
|                          |        | n_h/n_e  | 0.20   |     | E_c/T_e   | 14.8000 | lambda_h | 0.2009   |
| T_e/KeV                  | 4.00   |          |        |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q        | 2.00   |     | tau_c     | 0.3482  |          |          |
|                          |        |          |        |     |           |         | nrx      | 174.4104 |
| E_0/KeV                  | 170.00 | e_i      | 1.00   |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h      | 1.00   |     |           |         |          |          |
| a                        | 1.26   |          |        |     | n_i/n_e   | 0.8000  |          |          |
| R                        | 6.29   |          |        |     |           |         |          |          |
| B0/T                     | 1.3750 |          |        |     |           |         |          |          |
|                          |        |          |        |     |           |         |          |          |
| perturbation wave length | 0.90   |          |        |     | 参量（程序调参）  |         |          |          |
|                          |        |          |        |     | T_h/KeV   | 8500    |          |          |
|                          |        |          |        |     | v_0/C_s^h | 0.2     |          |          |
0.04-0.96 DeltaLambdasdf = 0.05

| cr | ω1       | γ1       | ω2       | γ2       |
|----|----------|----------|----------|----------|
| 理论 | 1.25738  | 0.15068  | 3.66852  | 0.00883  |
| 模拟 | 1.25019  | 0.13359  | 3.57531  | 0.01979  |
| 误差 | 0.57%    | 11.34%   | 2.54%    | 124.20%  |

### t262 13 NLT_1 2146326
[[#t260 13 2143658|t260]]
time * 3
## 2024-03-25
### t263 13 NLT 2162551
[[#t261 13 NLT 2144228|t261]]
n_h=0.55
n_i=0.45

| 输入参数                     |         |          |         |  |           |          |          |           |
|--------------------------|---------|----------|---------|--|-----------|----------|----------|-----------|
| m_h                      | 1.00    | Lambda 0 | 0.3000  |  | 参量（理论对比）  |          |          |           |
| m_i                      | 1.00    |          |         |  | E_0/T_e   | 42.5000  | lambda_i | 0.0368    |
|                          |         | n_h/n_e  | 0.55    |  | E_c/T_e   | 14.8000  | lambda_h | 0.2009    |
| T_e/KeV                  | 4.00    |          |         |  |           |          |          |           |
| T_i/KeV                  | 4.00    | q        | 2.00    |  | tau_c     | 0.3482   |          |           |
|                          |         |          |         |  |           |          | nrx      | 174.4104  |
| E_0/KeV                  | 170.00  | e_i      | 1.00    |  | E_c/KeV   | 59.2000  |          |           |
|                          |         | e_h      | 1.00    |  |           |          |          |           |
| a                        | 1.26    |          |         |  | n_i/n_e   | 0.4500   |          |           |
| R                        | 6.29    |          |         |  |           |          |          |           |
| B0/T                     | 1.3750  |          |         |  |           |          |          |           |
|                          |         |          |         |  |           |          |          |           |
| perturbation wave length | 0.90    |          |         |  | 参量（程序调参）  |          |          |           |
|                          |         |          |         |  | T_h/KeV   | 8500     |          |           |
|                          |         |          |         |  | v_0/C_s^h | 0.2      |

| cr | ω1       | γ1       | ω2       | γ2       |
|----|----------|----------|----------|----------|
| 理论 | 0.87374  | 0.08839  | 4.49729  | 0.04121  |
| 模拟 | 0.87065  | 0.08479  | 4.30300  | 0.05228  |
| 误差 | 0.35%    | 4.07%    | 4.32%    | 26.86%   |

### t264 13 NLT 2164549
[[#t261 13 NLT 2144228|t261]]
[t264](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct264)

| 输入参数                     |        |             |          |     |           |         |          |          |
| ------------------------ | ------ | ----------- | -------- | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0    | 0.3000   |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |             |          |     | E_0/T_e   | 42.5000 | lambda_i | 0.0368   |
|                          |        | ==n_h/n_e== | ==0.65== |     | E_c/T_e   | 14.8000 | lambda_h | 0.2009   |
| T_e/KeV                  | 4.00   |             |          |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q           | 2.00     |     | tau_c     | 0.3482  |          |          |
|                          |        |             |          |     |           |         | nrx      | 174.4104 |
| E_0/KeV                  | 170.00 | e_i         | 1.00     |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h         | 1.00     |     |           |         |          |          |
| a                        | 1.26   |             |          |     | n_i/n_e   | 0.3500  |          |          |
| R                        | 6.29   |             |          |     |           |         |          |          |
| B0/T                     | 1.3750 |             |          |     |           |         |          |          |
|                          |        |             |          |     |           |         |          |          |
| perturbation wave length | 0.90   |             |          |     | 参量（程序调参）  |         |          |          |
|                          |        |             |          |     | T_h/KeV   | 8500    |          |          |
|                          |        |             |          |     | v_0/C_s^h | 0.2     |          |          |

| cr  | ω1      | γ1      | ω2      | γ2      |
| --- | ------- | ------- | ------- | ------- |
| 理论  | 0.78881 | 0.04841 | 4.65396 | 0.04785 |
| 模拟  | 0.78393 | 0.05626 | 4.43608 | 0.05692 |
| 误差  | 0.62%   | 16.22%  | 4.68%   | 18.97%  |

### t265 13 NLT_1 2166625

| 输入参数                     |        |             |          |     |           |         |          |          |
| ------------------------ | ------ | ----------- | -------- | --- | --------- | ------- | -------- | -------- |
| m_h                      | 1.00   | Lambda 0    | 0.3000   |     | 参量（理论对比）  |         |          |          |
| m_i                      | 1.00   |             |          |     | E_0/T_e   | 42.5000 | lambda_i | 0.0366   |
|                          |        | ==n_h/n_e== | ==0.65== |     | E_c/T_e   | 14.8000 | lambda_h | 0.1997   |
| T_e/KeV                  | 4.00   |             |          |     |           |         |          |          |
| T_i/KeV                  | 4.00   | q           | 2.00     |     | tau_c     | 0.3482  |          |          |
|                          |        |             |          |     |           |         | nrx      | 175.8657 |
| E_0/KeV                  | 170.00 | e_i         | 1.00     |     | E_c/KeV   | 59.2000 |          |          |
|                          |        | e_h         | 1.00     |     |           |         |          |          |
| a                        | 0.98   |             |          |     | n_i/n_e   | 0.3500  |          |          |
| R                        | 3.75   |             |          |     |           |         |          |          |
| B0/T                     | 2.0000 |             |          |     |           |         |          |          |
|                          |        |             |          |     |           |         |          |          |
| perturbation wave length | 0.80   |             |          |     | 参量（程序调参）  |         |          |          |
|                          |        |             |          |     | T_h/KeV   | 8500    |          |          |
|                          |        |             |          |     | v_0/C_s^h | 0.2     |          |          |
delta lambda=0.05               0.04-0.86 

| cr  | ω1      | γ1      | ω2      | γ2      |
| --- | ------- | ------- | ------- | ------- |
| 理论  | 0.78879 | 0.04834 | 4.65653 | 0.04745 |
| 模拟  | 0.78455 | 0.05295 | 4.42964 | 0.05481 |
| 误差  | 0.54%   | 9.53%   | 4.87%   | 15.50%  |
**符合**理论#

### t266 13 NLT 2166723
[[#t264 13 NLT 2164549|t264]]
time * 3


### t267 13 NLT_2 2167495
[[#t251 13 170 3432973 chenzhe|t251]]
delta lambda=0.05
**修改EP密度分布** 居中 [[lu2019theoretical]]
对比电场变化

| cr  | ω1      | γ1      | ω2      | γ2      |
| --- | ------- | ------- | ------- | ------- |
| 理论  | 1.25732 | 0.15069 | 3.66918 | 0.00873 |
| 模拟  | 1.25502 | 0.14289 | 3.49858 | 0.01801 |
| 误差  | 0.18%   | 5.17%   | 4.65%   | 106.23% |
可能因为γ2太小导致误差不可避免很大

## 2024-03-26
### t268 13 NLT_2 2173376
[[#t265 13 NLT_1 2166625|t265]]
**修改EP密度分布** 居中 rwf=0.02

| cr  | ω1      | γ1      | ω2      | γ2      |
| --- | ------- | ------- | ------- | ------- |
| 理论  | 0.78879 | 0.04834 | 4.65653 | 0.04745 |
| 模拟  | 0.72863 | 0.04421 | 4.11775 | 0.13579 |
| 误差  | 7.63%   | 8.55%   | 11.57%  | 186.16% |

### t269 13 NLT_2 2173405
[[#t265 13 NLT_1 2166625|t265]]
```
**修改EP密度分布**  rwf=0.0002
```

| cr  | ω1      | γ1      | ω2      | γ2      |
| --- | ------- | ------- | ------- | ------- |
| 理论  | 0.78879 | 0.04834 | 4.65653 | 0.04745 |
| 模拟  | 0.96469 | 0.08145 | 4.13420 | 0.08573 |
| 误差  | 22.30%  | 68.50%  | 11.22%  | 80.66%  |

### t270 13 NLT_3 2173944
减小初始扰动密度大小[[#t251 13 170 3432973 chenzhe|t251]]

| cr  | ω1      | γ1       | ω2      | γ2         |
| --- | ------- | -------- | ------- | ---------- |
| 理论  | 1.25732 | 0.150689 | 3.66918 | 0.00873177 |
| 模拟  | 1.25260 | 0.13586  | 3.61997 | 0.02990    |
| 误差  | 0.38%   | 9.84%    | 1.34%   | 242.42%    |

## 2024-03-27
### t271 13 NLT_3 2180336
Nh=0.65
Ni=0.35
[[#t270 13 NLT_3 2173944|t270]]

| cr  | ω1      | γ1      | ω2      | γ2      |
| --- | ------- | ------- | ------- | ------- |
| 理论  | 0.78879 | 0.04834 | 4.65653 | 0.04745 |
| 模拟  | 0.78798 | 0.05291 | 4.44628 | 0.04905 |
| 误差  | 0.10%   | 9.46%   | 4.52%   | 3.37%   |

### t272 13 NLT_3 2180401
```
**修改EP密度分布** rwf=0.0002
减小初始扰动大小
```

| cr | ω1       | γ1       | ω2       | γ2       |
|----|----------|----------|----------|----------|
| 理论 | 0.78879  | 0.04834  | 4.65653  | 0.04745  |
| 模拟 | 0.96575  | 0.08430  | 4.14097  | 0.10263  |
| 误差 | 22.43%   | 74.38%   | 11.07%   | 116.30%  |

## 2024-03-29
### t273 13 NLT_4 2203160
温度各向异性 

### t274 13 NLT_4 2203717
```
温度各向异性 ANTFM
**修改EP密度分布** rwf=0.0002
减小初始扰动大小
```

$2 \pi/\lambda=k$
$\lambda \sim 0.05a$
## 2024-03-30
### t275 13 NLT_2 2218676
[[#t265 13 NLT_1 2166625|t265]]
```
**修改EP密度分布**  rwf=0.0002
random init
```
## 2024-04-04

### t276 13 NLT_2 2281737
[[#t267 13 NLT_2 2167495|t267]] 
```markdown
**修改EP密度分布**  rwf=0.0002
random init
nlnr=T
maxstep=100000
simtitaltime=40.2
nalpha=2 -> nalpha = 382
```
[t276](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct276)
## 2024-04-06
### t277 13 NLT_2 2281872
```markdown
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=100000
simtitaltime=40.2
nalpha=2
```
[t277](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct277)
GAM理论652 -0.00607605 I
图中为1.5

### t278 13 NLT_2 2281905
```markdown
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=100000
simtitaltime=40.2
nalpha=2
q=2+0.4(r/a)^2
```
[t278](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct278)
q增大后，GAM位置更加趋向于芯部，宽度变窄
**GAM频率理论为1.8，图中为1.1左右**
电场幅值变大
## 2024-04-07
### t279 13 NLT_2 2292740
```markdown
**修改EP密度分布**  rwf=0.0002  rcf=0.5
**热离子**密度  与EP分布相反x
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2
```

### t280 13 NLT_2 2292868
```markdown
**修改EP密度分布**  rwf=0.0002  rcf=0.5
**热离子**密度  与EP分布相反x
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2+0.4(r/a)^2
```

### t281 13 NLT_2 2292881
```markdown
**修改EP密度分布**  rwf=0.0002  rcf=0.5
**热离子**密度  与EP分布相反x
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=3+0.4(r/a)^2
```


### t282 13 NLT_2 2292922
```markdown
**修改EP密度分布**  rwf=0.0002  rcf=0.5
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=3+0.4(r/a)^2
```

### t283 13 NLT_2 2296721

```
**修改EP密度分布**  rwf=0.0002  rcf=0.5
**热离子**密度  与EP分布相反
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=3+0.4(r/a)^2
```
[t283](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct283)
无高能离子时GAM:1.72595 -0.00217656 I
EGAM：0.512926 + 0.08638 I
3.16564 + 0.0542216 I
EGAM在频谱图上显示不清晰，低频分支消失，出现一支频率为2.2的波动
### t284 13 NLT_2 2296812
[[#t278 13 NLT_2 2281905|t278]]
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=3+0.4(r/a)^2
```
以一支频率为2.4的波动为主，其余不明显，不是GAM，不是EGAM
理论GAM 1.72194 -0.00178353 I
连续谱?
### t285 13 NLT_2 
高能离子密度扫描
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
**q=3+0.4(r/a)^2**
```
[testslurm.xlsx](file:///F:%5CSimulation%5Ctestslurm.xlsx)
[t285](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct285)

**nh=0.04** 0.8
**0.3a**GAM增长
## 2024-04-08
### t286 13 NLT_2 2310975
```
**修改EP密度分布**  rwf=0.1  rcf=0.2
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=3+0.4(r/a)^2
```
ni=0.7
nh=0.3


### t287 13 NLT_2 2311001
```
**修改EP密度分布**  rwf=0.1  rcf=0.2
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2
```
ni=0.7
nh=0.3
**flat q 干扰信号减少，高能离子浓度梯度引起连续谱**（local？）

### t288 13 NLT_2 2314291
```
**修改EP密度分布**  rwf=0.1  rcf=0.2
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2
```
ni=0.35
nh=0.65
**0.5处GAM信号为增长**
## 2024-04-09
### t289 13 NLT_NL 2319098
[[#t288 13 NLT_2 2314291|t288]] 
```
**修改EP密度分布**  rwf=0.1  rcf=0.2
random init
nlnr=T
maxstep=100000
simtitaltime=40.2
nalpha=286
q=2
nmu=96
```

```bash
error! nFltr < npro
 please set mod(nFltr,npro) = 0
 nFltr = 1 + FLOOR((nalpha-1.0)/nkc) =           1
 nalpha =           2
 nkc =           3
 npro =          96
```

非线性结果如何解释

linear:
0.99  
1  2319945
2 2320179



### t290 13 NLT_2 2319961
```
**修改EP密度分布**  rwf=0.1  rcf=0.1
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2
```
ni=0.35
nh=0.65


### t291 13 NLT_2 2319967
```
**修改EP密度分布**  rwf=0.1  rcf=0.3
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2
```
ni=0.35
nh=0.65


### t292 13 NLT_2 2320405
```
**修改EP密度分布**  rwf=0.1  rcf=0.6
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2
```
ni=0.35
nh=0.65

### t293 13 NLT_2 2320513
```
**修改EP密度分布**  rwf=0.1  rcf=0.6
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2
nrx=192--> 240
0.82->0.98
```
ni=0.35
nh=0.65


### t294 13 NLT_2 2320572
```
**修改EP密度分布**  rwf=0.1  rcf=0.3
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2
nrx=192--> 240
0.82->0.98
```
ni=0.35
nh=0.65

## 2024-04-11
### t295 13 NLT_NL 2342808
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=T
maxstep=100000
simtitaltime=40.2
nalpha=286
q=2
nmu=96
```
范围0.98
nh0.2  ni0.8
nrx=220
### t296 13 NLT_2 2339250
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2
nrx=192--> 240
```
nh 0.2 ni 0.8
范围0.98

### t297 13 NLT_2 2339283
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2
nrx=192--> 240
```
nh 0.2 ni 0.8
范围0.98
ion PType=1

### t298 13 NLT_2 2339288
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2
nrx=192--> 240
```
nh 0.2 ni 0.8
范围0.98
ion PType=2


## 2024-04-12
### t299 13 NLT_5 2359037
[[#t285 13 NLT_2|t285]]
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
**GAM init**
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=3+0.4(r/a)^2
```

n_h=0
GAM 初始化为阻尼
随机扰动且n_h为零时有增长的波动
### t300 13 NLT_6 2359064
[[#t298 13 NLT_2 2339288|t298]]
对比[[#t296 13 NLT_2 2339250|296]]q梯度
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2+0.5(r/a)
nrx=192--> 240
```
nh 0.2 ni 0.8
范围0.98

### t301 13 NLT_2 2359564
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=100000
simtitaltime=40.2
nalpha=2
q=2
nrx=192--> 240
```
nh 0.2 ni 0.8
范围0.98

### t302 13 NLT_2 2359567
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=100000
simtitaltime=40.2
nalpha=2
q=2
nrx=192--> 240
```
nh 0.2 ni 0.8
范围0.98
ion PType = 2

### t303 13 NLT_2 2493603
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=100000
simtitaltime=40.2
nalpha=2
q=1.75
nrx=192--> 240
```
nh 0.2 ni 0.8
范围0.98

### t304 13 NLT_2 2507377
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
GAM init
nlnr=F
maxstep=100000
simtitaltime=40.2
nalpha=2
q=2
nrx=192--> 240
RH扰动0.8-0.9
模拟范围0.01-0.99
```
nh = 0.2

nh = 0.3  2511685
### t305 NLT 对比GAM理论模拟
[[#test54 modi ed|test54]]
```
q=1  2509025
q=5  2509023
```

```
k = 0.203
q =1 2513054
q =5 2513061
```

### t306 mesh_F0 增长GAM NLT_1
[[#t288 13 NLT_2 2314291|t288]]
1 画F_0   2516501
```
**修改EP密度分布**  rwf=0.1  rcf=0.2
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=2
ni=0.35
nh=0.65
```
2 增加EP宽度 2516507
rcf = 0.3

3 改变浓度 2516512
ni=0.7
nh=0.3
rcf=0.2

## 2024-05-06
### t307 NLT_2 2524311
[[#t285 13 NLT_2]]
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=3
```

### t308 NLT_2 2534548
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=3.2564
```

### t309 NLT_2 2542359
```
**修改EP密度分布**  rwf=0.1  rcf=0.5
random init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
**q=3+0.4(r/a)^2**
```
nh=0


### t310 NLT_1(global GAM) 8673849

```
GAM init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
**q=3
```
nh=0 ni=1
R_L=0.85 R_R=0.90
0.1-0.99

### t311 NLT_2
```
RANDOM init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=3+0.4*(r/a)**2
```
nh=0 ni=1
slowing down[87548_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct311%5C87548_GAM_KE)
bump on tail [87550_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct311%5C87550_GAM_KE)
数据有区别，原因？
 重复
 slowing down 8902413
bump on tail 8902423
均为弱阻尼或阻尼
### t312 NLT_2
```
RANDOM init
nlnr=F
maxstep=50000
simtitaltime=20.1
nalpha=2
q=3+0.4*(r/a)**2
```
nh=0.04 ni=0.96
slowing down [96381_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct312%5C96381_GAM_KE)
bump on tail [96466_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct312%5C96466_GAM_KE)

## 2024-05-28
### t313 NLT_2
非线性
8907171
8907574
8912401
8912583 √
## tr NLT_rotation 0.1 
[tr](file:///F:%5CSimulation%5Crotationbenchmark%5Ctr)
### 1

```
R = 3.75
a = 0.4
B0 = 2.0
q = 2
L_R = 0.8
R_R = 0.9
T_i = 0.12
T_e = 0.0001
k_r rho_i = 0.1036
nrx = 80
```

### 2 2505162
```
R = 3.75
a = 0.4
B0 = 2.0
q = 2
L_R = 0.8
R_R = 0.9
T_i = 0.12
T_e = 0.0001
k_r rho_i = 0.1036
r_max = 0.99
r_min = 0.01
M = 0
nrx = 528
PType=4
```
模拟范围可以远超过扰动的长度，观察大径向尺度内GAM的演化
output：除扰动范围内有GAM，径向接近0的位置的幅度增长如何解释？

### 3 2508062
```
R = 3.75
a = 0.4
B0 = 2.0
q = 2
L_R = 0.8
R_R = 0.9
T_i = 0.12
T_e = 0.0001
k_r rho_i = 0.1036
r_max = 0.99
r_min = 0.01
M = -0.4
nrx = 528
PType=7
```
**为什么计算得到的初始电场值有区别**


## 2024-05-30
### t314 NLT_1
[[ren2020global]]
Global GAM
[29020_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct314%5C29020_GAM_KE)重复：[54805_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct314%5C54805_GAM_KE)  [68455_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct314%5C68455_GAM_KE)
改变扰动位置

[30662_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct314%5C30662_GAM_KE) [26493_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct314%5C26493_GAM_KE)


[54817_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct314%5C54817_GAM_KE) 

观察到GAM连续谱及不随径向位置改变的频率模式，两个模式频率谱强度均在扰动位置处最大，
“全局”模式在扰动位置处与连续谱交点出现强度峰值，但不是连续谱频率最大值点

## 2024-06-01
### t315 NLT_1
不加温度梯度 [68486_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct315%5C68486_GAM_KE)


## 2024-06-02
### t316 NLT_1
具有离轴最大值的连续谱
[74259_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct316%5C74259_GAM_KE)
结果不明显
### t317 NLT_1
随机扰动，去除RHtest  [74857_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct317%5C74857_GAM_KE)
仅出现连续谱，峰值明显离轴
### t318 NLT_1
[75154_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct318%5C75154_GAM_KE)
具有离轴最大值的连续谱
RHtest

## 2024-06-04
### t319 NLT_1
flat profile
[97726_GAM_KE](file:///media/imyxl/T7%20Shield/Simulation/biaozhun/EPs1.0/t319/97726_GAM_KE)


## 2024-06-05
### t320 NLT_1
尝试消除边界数据不稳定 
[41770_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct320%5C41770_GAM_KE)

## 2024-06-09
### t321 NLT_1
[03773_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct321%5C03773_GAM_KE)

## 2024-06-20
### t322 NLT_1
连续谱无极值点[89845_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct322%5C89845_GAM_KE)

### t323 NLT_1
连续谱无极值点[91448_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct323%5C91448_GAM_KE)

## 2024-06-21
### t324 NLT_1
[[miyato2006nonlocal]]
GAM_initial  [57009_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct324%5C57009_GAM_KE)
random_initial [57059_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct324%5C57059_GAM_KE)

## 2024-07-04
### t325 NLT_1
[[kong2013observation]]   [[li2011investigation]]

random_initial[00992_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct325%5C00992_GAM_KE)
密度减小 [01241_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct325%5C01241_GAM_KE)


# 本征模激发
## 2024-08-18
### t326 NLT_1
9826772 \9829525 测试本征值激发频率宽度
9829866  四个参考点正弦扰动
9830086 八个参考点
模拟解大致范围:

数值解大致范围?

如何获得模拟解模结构？

## 2024-08-19
### t327 NLT_1
(1)
GAM_B init + random_init 9830821[30821_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct327%5C30821_GAM_KE)
random_init 9830830[30830_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct327%5C30830_GAM_KE)

(2)
GAM_B init + random_init 9844635 [44635_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct327%5C44635_GAM_KE)
random_init 9844628 [44628_GAM_KE](file:///F:%5CSimulation%5Cbiaozhun%5CEPs1.0%5Ct327%5C44628_GAM_KE)
对比结果：激发后二阶向外传
### t328 NLT_1
采用数值解剖面

### t329 NLT_1
[[xiang2020weilong]]（EGAM）
验证径向电场是否为单一频率的全局模结构，还是连续谱组成结构
1.GAM模拟 10550475
总结不出模结构特点

2.GAM tau =1 10550693
内部模结构为平的，中间向外为振荡结构

3.EGAM 10553917
原GAM的模结构结果之外出现了文章中EGAM典型的模结构，在内部出现了振荡，说明该结果仅由高能粒子引起



### 25_01_01_Lambda_scan

11206336
delta lambda = 0.2

11206343
delta lambda = 0.5

11206346
delta lambda = 1

11206349
delta lambda = 2

11206357
delta lambda = 5