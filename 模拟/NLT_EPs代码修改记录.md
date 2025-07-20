## eq_adhoc.f90
```fortran
!!! eq_adhoc.f90
**flat q profile**
!304
q_bar(i) = (c0 + c1*(r(i)/a) + c2*(r(i)/a)**2 + c3*(r(i)/a)**3 + c4*(r(i)/a)**4)*sqrt(1._dp-(r(i)/R0)**2)
!411
qbt = (c0 + c1*(tmpr/a) + c2*(tmpr/a)**2 + c3*(tmpr/a)**3 + c4*(tmpr/a)**4)*sqrt(1._dp-(tmpr/R0)**2)
```


## Filter.f90
```fortran
!!!Filter.f90
!198 q死循环

!$OMP PARALLEL DO PRIVATE(i)

    do i = 1, nrx

      !n = (nFltr-1) * ntm

      !tmp = 0._dp

      !print*,nFltr,ntm,dmFltr

      !do while (tmp .lt. mnFactorC)

       ! n = n - ntm

       ! do m = NINT(n*mesh_q(i)-dmFltr), NINT(n*mesh_q(i)+dmFltr)

        !  call comput_filter_factor(m,n,i,tmp)

         ! if (tmp .ge. mnFactorC) exit

        !end do

      !end do

      !n = n/ntm

      h_Lambda_alpha_c(i) = 10._dp

    end do

    !$OMP END PARALLEL DO
```

## Global.f90
### ed1
```fortran
!!!Global.f90

!slowing down EPs
  real(dp) :: SDFfactor,Lambda0sdf,DeltaLambdasdf,critienerfasdf
  integer :: enerpariPar,redefinedmt

!perturbation
  real(dp) :: r_max,r_min,lalph,R_R,L_R

  namelist /PARTICLE/  m_ref_new,IonNum,FlagKE,N_ref,T_ref,Lambda0sdf,DeltaLambdasdf

  namelist /GRID/      FlagAxis,r_max,r_min,rxType,nrx,ntm,nalpha,ntheta, &
nvpara,muType,nmu,dtMax,simTotalTime,maxStep,R_R,L_R

!归一化质量、温度；
    redefinedmt = 1
    call ReadParticle
    m_ref_new = Parlist(1)%m*m_ref_new
    T_ref = Parlist(1)%TPar*T_ref
    redefinedmt = 2

!slowing down distribution function factor
    do iPar = 1, IonNum
      if (ParList(iPar)%PT .eq. 5 .or. ParList(iPar)%PT .eq. 6) then
        enerpariPar = iPar
      end if
    end do
    if (Parlist(enerpariPar)%PT .eq. 6) then
      do iPar = 1, IonNum
        critienerfasdf = critienerfasdf+Parlist(iPar)%NPar*(Parlist(iPar)%e)**2/(Parlist(iPar)%m*m_ref_new)
      end do
      critienerfasdf = (critienerfasdf/Parlist(0)%NPar)**(2._dp/3._dp)*Parlist(enerpariPar)%m*m_ref_new
      SDFfactor = integrate_sdf(Lambda0sdf,DeltaLambdasdf)
    end if

!0.9->1.0
if (r_max .ge. 1.0) then

!归一化质量、温度
    if (redefinedmt .eq. 1) then
    allocate(ParList(0:IonNum))
    end if

!Double shifted maxwell distribution function
  !energetic particles PType=5
  function FDSHM(i,k,l,ms,mu,N,U,T)
    implicit none
    integer :: i,k,l
    real(dp) :: ms,v1,v2,mu,B,N,T,U,FDSHM
    v1 = Parlist(iPar)%Vvpara0(l) - U
    v2 = Parlist(iPar)%Vvpara0(l) + U
    B = mesh_B(i,k)
    FDSHM = 0.5_dp*N*(sqrt(ms/(twopi*T))**3)*exp(-(mu*B)/T)*(exp(-(0.5_dp*ms*v1*v1)/T)+exp(-(0.5_dp*ms*v2*v2)/T))
  end function FDSHM
  !energetic particles

!slowing down distribution and integrate
    !energetic particles PType=6
  function SDF(i,k,l,ms,mu,N,U)
    implicit none
    integer :: i,k,l
    real(dp) :: ms,mu,B,N,U,SDF,vtotal,kinetener,Lambda,c,vc
    B = mesh_B(i,k)
    vtotal = sqrt((Parlist(iPar)%Vvpara0(l))**2+2._dp*mu*B/ms)
    kinetener = 0.5_dp*ms*Parlist(iPar)%Vvpara0(l)**2+mu*B
    Lambda = mu*B0/kinetener
    if (U.ge.vtotal) then
      vc = sqrt(2.0_dp*14.8_dp*Parlist(0)%Tt(i)*critienerfasdf/ms)
      c = 2.0_dp*PI*(log(U**3+vc**3)-log(vc**3))*SDFfactor/3.0_dp
      c = 1.0_dp/c
      SDF = c*N*exp(-((Lambda-Lambda0sdf)/DeltaLambdasdf)**2)/(vtotal**3+(vc)**3)
    else
      SDF = 0.0_dp
    end if
  end function SDF
  !energetic particles

  function integrate_sdf(la0, dla)
    real(dp), intent(in) :: la0, dla
    real(dp) :: la, h,integral,integrate_sdf
    integer :: i, n

    n = 10000   ! 设置积分分割的数量
    h = 1.0_dp / real(n)
    integral = 0.0_dp

    !$OMP PARALLEL DO PRIVATE(i, la) reduction(+:integral)
    do i = 1, n-1
      la = real(i) * h
      integral = integral + exp(-((la - la0)/dla)**2) / sqrt(1.0_dp - la)
    end do
    !$OMP END PARALLEL DO

    integral = h * (0.5_dp * (exp(-(la0/dla)**2) / sqrt(1.0_dp - la0) + exp(-((1.0_dp - la0)/dla)**2) / sqrt(la0)) + integral + 0.5_dp * exp(-((1.0_dp - la0)/dla)**2) / sqrt(1.0_dp - la0))

    integrate_sdf = integral
  end function integrate_sdf


!density,EPs
    else if (PType .eq. 6) then
      rm     = 0.40_dp
      deltaN = 0.05_dp
      kappaN = 2.2320_dp
      DenN0 = exp(-kappaN*ra/R0*( r - rm &
                   - deltaN*tanh((r-Vmr(1)  )/deltaN) &
                   - deltaN*tanh((r-Vmr(nrx))/deltaN) ))**对于高能粒子分布不合理**


```


### ed2
```Fortran
real(dp) :: cnf1,rcf,rwf

    if (PType .eq. 1) then
      rm     = 0.50_dp
      deltaN = 0.30_dp
      kappaN = 2.23_dp
      DenN0 = exp(-kappaN*ra/R0*deltaN*tanh((r-rm)/deltaN))
    else if (PType .eq. 2) then
      rm     = 0.40_dp
      deltaN = 0.05_dp
      kappaN = 2.2320_dp
      DenN0 = exp(-kappaN*ra/R0*( r - rm &
                   - deltaN*tanh((r-Vmr(1)  )/deltaN) &
                   - deltaN*tanh((r-Vmr(nrx))/deltaN) ))
    else if (PType .eq. 6) then
      cnf1=0.5_dp
      rcf=0.4_dp
      rwf=0.02_dp
      DenN0=(1.0_dp+cnf1*(tanh(rcf/rwf)-1.0_dp))**(-1.0_dp)*(1.0_dp+cnf1*(tanh((rcf-r)/rwf)-1.0_dp))
```
## mesh.dat
```fortran
!!mesh.dat

      Lambda0sdf = 0.50
    DeltaLambdasdf = 0.05

!
           R_R = 0.60 
           L_R = 0.45 
```

## VlasovDF.f90
```fortran
!!!VlasovDF.f90

    !real(dp),parameter :: L_R = 0.45_dp, R_R = 0.60_dp

!Random Init
          if (.FALSE.) then


!GAM_B Init
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

```

## VlasovF0.f90
```fortran
!!!VlasovF0.f90

!输出速度网格，分布函数
      do i = 0, nmu - 1
        if (pro_id .eq. i) then
          ! 构造文件名，包括 iPar 和 pro_id
          write(unit=filenamedata, fmt="('grid2dF0_',I2.2,'_pro',I2.2,'.dat')") iPar, pro_id
          filenamedata = trim(filenamedata)
          
          ! 使用关联的文件单元号打开文件
          open(111, file=filenamedata, status='unknown', position='append')
          
          do l = 1, nvpara
            write(111, *) ParList(iPar)%Vvpara0(l), ParList(iPar)%Vmu0(i+1), ParList(iPar)%F0(i_diag,k_diag,l)
          end do
          
          ! 关闭文件
          close(111)
        end if
      end do

!energetic particles
    if (iPar .ne. enerpariPar) then
      !$OMP PARALLEL DO PRIVATE(i,k,l)
      do l = 1, nvpara
      do k = 1, ntheta
      do i = 1, nrx
        ParList(iPar)%F0(i,k,l) = FM(i,k,l,m_s,mu_loc,ParList(iPar)%Nt(i), &
                                  ParList(iPar)%Ut(i),ParList(iPar)%Tt(i))
      end do
      end do
      end do
      !$OMP END PARALLEL DO

     !energetic particles PType=5:double shifted ;PType=6:slowing down distribution
    else if (Parlist(iPar)%PT .eq. 5) then
      !$OMP PARALLEL DO PRIVATE(i,k,l)
      do l = 1, nvpara
      do k = 1, ntheta
      do i = 1, nrx
        ParList(iPar)%F0(i,k,l) = FDSHM(i,k,l,m_s,mu_loc,ParList(iPar)%Nt(i), &
                                  ParList(iPar)%Ut(i),ParList(iPar)%Tt(i))
      end do
      end do
      end do
      !$OMP END PARALLEL DO
    else if (Parlist(iPar)%PT .eq. 6) then
      !$OMP PARALLEL DO PRIVATE(i,k,l)
      do l = 1, nvpara
      do k = 1, ntheta
      do i = 1, nrx
        ParList(iPar)%F0(i,k,l) = SDF(i,k,l,m_s,mu_loc,ParList(iPar)%Nt(i), &
                                  ParList(iPar)%Ut(i))
      end do
      end do
      end do
        !$OMP END PARALLEL DO
    end if
    !energetic particles
```


## **未修改文件**
Equilibrium.f90
Diagnoses.f90
Field.f90
FieldA.f90
FieldB.f90
GyroAverage.f90
Initial.f90
InrType.f90
mklDfti.f90
NLT.f90
Numerical.f90
Source.f90
SourceA.f90
SourceB.f90
Spline.f90
Test.f90