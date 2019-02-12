module globals
  !一般変数 for general
 integer,save:: i,j,k,mcs,iseed,kt,ks
 integer,save:: index_n
 real(8),save:: T,beta
 character(80),save:: f2
  !パラメーター parameter
 integer,parameter :: L = 24
 integer,parameter :: N = L * L
 integer,parameter :: S = 10
 real(8),parameter :: pi = 2.0d0 * acos(0.0d0)
 integer,parameter :: NT = 200
 real(8),parameter :: dT = 0.025d0
 integer,parameter :: Nmcs = 30000
 integer,parameter :: Mcs_min = 1000
  integer,parameter:: Neff = Nmcs - Mcs_min
 real(8),parameter :: T_ini = 5.0d0
 character(80),save :: f3(NT),f1
  !スピン for_Spin_Configuration
 real(8),save::Sx(N),Sy(N)
  !メトロポリス for_Metroplis_sweep
 real(8),save:: hx,hy
 real(8),save:: ep,ea,delta
 real(8),save:: Sx_d(N),Sy_d(N)
 real(8),save:: nx,ny,tempx,tempy,temp2x,temp2y,temp
 real(8),save:: cos_alpha,alpha
  !物理量 for_physical_quantities
 real(8),save:: Magx,Magy,Mag_all,Mag,Mag2
 real(8),save:: de,Energy,Energy1,Energy2
 real(8),save:: E(NT,S),M(NT,S),C(NT,S),X(NT,S)
  !相互作用の相手 for_partner_of_interaction
 integer,save:: ii_1(N),ii_2(N),ii_3(N),ii_4(n)
 real(8),save:: position(L,L)
 real(8) theta(N)
  !周期境界条件 for_preiodic_condition
 integer,save:: p(0:L+1)
  !サンプル平均 for sample_average
 real(8),save:: E_ave,M_ave,C_ave,X_ave
 real(8),save:: sd_E,sd_M,sd_C,sd_X
 character(80),save:: n1,n2,n3,n4
end module globals
!
!**************************************************************************
!
module subprogs
  use globals
  use module_random
 implicit none
contains
 subroutine convert_from_2D_to_1D 
   index_n = 1
   do j = 1,L
      do i= 1,L
         position(i,j) = index_n
         index_n = index_n + 1
      enddo
   enddo
   index_n = 1
   do j = 1,L
      do i = 1,L
         ii_1(index_n) = position(p(i+1),j)
         ii_2(index_n) = position(p(i-1),j)
         ii_3(index_n) = position(i,p(j+1))
         ii_4(index_n) = position(i,p(j-1))
         index_n = index_n + 1
      enddo
   enddo
  end subroutine convert_from_2D_to_1D
!
!**********************************
!
subroutine name
  n1= 'Energy'
  n2= 'Magnetization'
  n3= 'Specific heat'
  n4= 'Susceptibility'
end subroutine name
!
!********************************　
!
 subroutine set_initial_spin_configuration
   call rfr(N)
   do i = 1,N
      theta(i) = 2.d0*pi*ur(i)
      Sx(i) = cos(theta(i))
      Sy(i) = sin(theta(i))
   enddo
  end subroutine set_initial_spin_configuration
!
!**********************************
!
 subroutine initialize_expectation_values
   de = 0.0d0
   Energy1 = 0.0d0
   Energy2 = 0.0d0
   Mag = 0.0d0
   Mag2 = 0.0d0
  end subroutine initialize_expectation_values
!
!************************************
!
 subroutine Metropolis
   nx = cos(2.0d0*pi*ur(mcs))
   ny = sin(2.0d0*pi*ur(mcs))
   do i = 1, N
      !分子場の計算 calculate_molecular_field
      hx = Sx(ii_1(i)) + Sx(ii_2(i)) + Sx(ii_3(i)) + Sx(ii_4(i))
      hy = Sy(ii_1(i)) + Sy(ii_2(i)) + Sy(ii_3(i)) + Sy(ii_4(i))
      ep = -1.0d0 * (Sx(i) * hx + Sy(i) * hy) !回転前の局所エネルギー
      cos_alpha = Sx(i) * nx + Sy(i) * ny
      alpha = acos(cos_alpha)
      Sx_d(i) = cos(2.0d0 * alpha) * Sx(i) - sin(2.0d0 * alpha) * Sy(i)
      Sy_d(i) = sin(2.0d0 * alpha) * Sx(i) + cos(2.0d0 * alpha) * Sy(i)
      ea = -1.0d0 * (hx * Sx_d(i) + hy * Sy_d(i)) !回転後の局所エネルギー
      if (ea-ep < 0)then
         Sx(i) = Sx_d(i)
         Sy(i) = Sy_d(i)
      else
         delta = -beta * (ea-ep)
         if (ur(i+Nmcs)< exp(delta)) then
            Sx(i) = Sx_d(i)
            Sy(i) = Sy_d(i)
         endif
      endif
   enddo
  end subroutine Metropolis
!
!**********************************
!
 subroutine physical_quantities_calculation_per_MCS
   do i = 1, N
      hx = Sx(ii_1(i)) +Sx(ii_2(i)) + Sx(ii_3(i)) + Sx(ii_4(i))
      hy = Sy(ii_1(i)) +Sy(ii_2(i)) + Sy(ii_3(i)) + Sy(ii_4(i))
      Energy = Energy - 0.50d0 * (hx * Sx(i) + hy * Sy(i))
   enddo
   Energy1 = Energy1 + Energy
   Energy2 = Energy2 + Energy * Energy
   Magx = sum(Sx(:))
   Magy = sum(Sy(:))
   Mag_all  = sqrt(Magx * Magx + Magy* Magy)
   Mag = Mag + Mag_all
   Mag2 = Mag2 + Mag_all * Mag_all
  end subroutine physical_quantities_calculation_per_MCS
!
!*******************************
!
 subroutine phsyical_quantities_after_montecarlo_sweep
   Energy1 = Energy1 / (N * Neff)
   Energy2 = Energy2 / (N * N * Neff)
   Mag = Mag / Neff
   Mag2 = Mag2 / Neff
   E(kt,ks) = Energy1
   M(kt,ks) = Mag
   C(kt,ks) = N * beta * beta * (Energy2 - Energy1 * Energy1)
   X(kt,ks) = beta * (Mag2 - Mag * Mag) / N
  end subroutine phsyical_quantities_after_montecarlo_sweep
!
!*******************************
! 
  subroutine sample_average(f2,T,T_ini,NT,S,A,A_ave,sd_A)
    character(80) f2
    real(8) T,T_ini
    integer kt,ks,NT,S
    real(8) A_ave,sd_A,A(NT,S)
    write(10,'(a\)') '##',f2
    write(10,*) ''
    T = T_ini
    do kt = 1,NT
       A_ave = 0.0d0
       sd_A  = 0.0d0
       do ks = 1,S
          A_ave = A_ave + A(kt,ks)
       enddo
       A_ave = A_ave / S
       do ks = 1,S
          sd_A = sd_A + (A(kt,ks) - A_ave) ** 2
       enddo
       sd_A = sqrt(sd_A / sqrt(dble(S*(S-1))))
       write(10,'(3e20.5)') T,A_ave,sd_A
       T = T - dT
    enddo
    write(10,*) ''
    write(10,*) ''
  end subroutine sample_average
!
!********************************
!
  subroutine Vector
    write(f3(kt),'("./Vector/","L",i4,"T",i4,".d")')1000+L,1000+kt
    open(kt,file=f3(kt))
    k = 1
    do j = 1,L
       do i = 1,L
          write(kt,'(2i10,2f15.3)') i,j,0.50d0 *Sx(k),0.50d0 * Sy(k)
          k=k+1
       enddo
    enddo
    close(kt)
  end subroutine Vector
end module subprogs
!
!****************************************************************
!
 
program main
  use globals
  use module_random
  use subprogs
  implicit none
  !周期境界条件の設定
  do i = 1,L
     p(i) = i
  enddo
  p(0)   = L
  p(L+1) = 1
 
  write(*,'(a\)') 'input iseed : '
  read(*,*) iseed
  call warmdrn(iseed,100)
  !配列の一次元化 convert_from_2D_to_1D
  call convert_from_2D_to_1D
  call name
  do ks = 1,S                   !サンプルループ sample_loop
     write(*,*) ks
     !初期のスピン配置 set_initial_spin_configuration
     call set_initial_spin_configuration
     T = T_ini
     do kt = 1,NT                !温度ループ temperture_loop
        beta = 1.0d0 / T
        call rfr(Nmcs+N)!乱数の内訳 1:Nmcs:回転ベクトル用,Nmcs+1:Nmcs+n:メトロポリス用
        call initialize_expectation_values
        do mcs = 1, Nmcs        !モンテカルロループ MC_loop
           Energy = 0.0d0
           !メトロポリス法
           call Metropolis
           if (mcs > Mcs_min) then
              call physical_quantities_calculation_per_MCS
           endif
        enddo
        if (ks == 1) then
           call Vector
        endif
        call phsyical_quantities_after_montecarlo_sweep
        T = T -dT
     enddo
     iseed = iseed + 1
  enddo
  write(f1,'("L",i4,"average.dat")')L
  open(10,file=f1)
  call sample_average(n1,T,T_ini,NT,S,E,E_ave,sd_E)
  call sample_average(n2,T,T_ini,NT,S,M,M_ave,sd_M)
  call sample_average(n3,T,T_ini,NT,S,C,C_ave,sd_C)
  call sample_average(n4,T,T_ini,NT,S,X,X_ave,sd_X)
  close(10)
end program main