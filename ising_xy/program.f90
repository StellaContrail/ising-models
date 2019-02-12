! define J = 1
module global
    implicit none
    double precision,parameter :: pi = acos(-1.0d0)
end module

module extension
    implicit none
    integer,private :: fdigit = 0
contains
    ! 出力処理(Any)が全て終わった後に呼び出す
    ! value / max : % (value out of max)
    ! isflushed : Flushする予定があるときに使う. 主にプログレスバーとともに出力結果も同時に出したいときに使う. 使わないときはfalseに設定する。
    subroutine pbout(value, max, isflushed)
        integer,intent(in) :: value, max
        double precision,save :: rate = 0d0, time = 0d0
        double precision dt, dr, estimate
        integer,parameter :: digit = 50
        integer remain(2), values(8), i
        character(len=20) :: FMT
        logical,intent(in) :: isflushed
        if (isflushed) then
            write (*, *)
        end if

        ! 変数maxの桁数を調べて書式指定子の作成を行う
        write (FMT, '(a, i0, a)') "(2(i", int(log10(real(max))) + 1, ",a))"

        ! 日時の取得
        call date_and_time(values=values)
        dr = rate   ! 前に実行したrateをひとまずdrに格納する
        dt = time   ! 前に実行したtimeをひとまずdtに格納する

        ! 時間の更新 : 今日の0時からどれだけ時間が経ったか. -> 一日を超える計算は正確に表示できない
        time = ((values(5)*60d0+values(6))*60d0+values(7))*1d3+values(8)
        dt = time - dt

        ! 割合の更新
        rate = dble(value) / dble(max)
        dr = rate - dr

        ! 残り時間の計算 (milliseconds)
        estimate = (1d0 - rate) * (dt / dr)

        ! min/sec表記に変換
        remain(1) = int(estimate/6d4) ! minutes
        remain(2) = int((estimate-remain(1)*6d4)*1d-3) ! seconds

        write (*, FMT, advance='no') value, " / ", max, " ["
        do i = 1, int(digit*rate)
            write (*, '(a)', advance='no') "="
        end do
        write (*, '(a)', advance='no') ">"
        do i = int(digit*rate) + 1, digit
            write (*, '(a)', advance='no') "-"
        end do

        if (isflushed) then
            fdigit = 2 * int(log10(real(max))) + 75
        end if
        write (*, '(a, f7.2, a, i3, a, i2, a, a)', advance='no') "]", 100d0*rate, "% ", remain(1), "m", remain(2), "s", char(13)

        if (value == max) then
            write (*, *)
        end if
    end subroutine

    ! 最初の出力処理(Any)が始まる前に呼び出す
    ! pbout()でisflushedをtrueにしたときのみ必要になる。falseのときはいらない。
    subroutine pbflush()
        character(len=9) FMT
        if (fdigit == 0) then
            return
        end if

        write (FMT, '(a, i0, a)') "(", fdigit, "x, 3a)"
        write (*, FMT, advance='no') char(13), char(8), char(13)
    end subroutine
end module

module subprog
    use global
    implicit none
contains
    ! inits : スピンの初期化を行う。すべてのスピンを+-1か0に設定する。
    ! init  : 0 or +-1
    ! s     : 各格子点のスピンの向きを格納する配列
    ! n     : 全格子点の数

    ! init = 0 -> All spins are set to random number
    ! init != 0 -> All spins are set to 0
    subroutine inits(init, s, n)
        integer, intent(in) :: n, init
        double precision, intent(out) :: s(n)
        integer i
        double precision rnd(n)
        call random_number(rnd)
        if (init == 0) then
            do i = 1, n
                s(i) = 2.0d0 * pi * rnd(i)
            end do
        else
            do i = 1, n
                s(i) = 0.0d0
            end do
        end if
    end subroutine inits
    
    ! ising : スピンの更新を行う。
    ! beta  : 1/Tを表す
    ! s     : すべての格子点のスピンの情報を格納する配列
    ! n     : すべての格子点の数
    ! hit   : 反転させたスピンの数
    ! ne    : それぞれの格子点の隣接する配位の番地を格納する配列
    subroutine ising(beta, s, n, hit, ne)
        double precision, intent(in) :: beta
        integer, intent(in) :: n, ne(n,4)
        double precision, intent(inout) :: s(n)
        double precision, intent(out) :: hit
        integer i, m, nhit
        double precision rnd1(n), rnd2(n), rnd3(n), old_energy, new_energy, s_f, ds
        
        nhit = 0
        ! 配列rnd1とrnd2に乱数を設定する
        call random_number(rnd1)
        call random_number(rnd2)
        call random_number(rnd3)

        do i = 1, n
            ! ランダムな格子点を選択し、番地番号に変換する
            m = int(dble(n)*rnd1(i)) + 1
            ! 稀に番地が範囲を超えることがあるのでその時はmに設定する
            if(m > n) then
                m = n
            endif
            ! ランダムな格子点をもし反転させたときのエネルギーを計算する。
            old_energy = -1.0d0*(cos(s(m)-s(ne(m, 1)))+cos(s(m)-s(ne(m, 2)))+cos(s(m)-s(ne(m, 3)))+cos(s(m)-s(ne(m, 4))))
            s_f = flip(s(m), 2d0*pi*rnd3(i))
            new_energy = -1.0d0*(cos(s_f-s(ne(m, 1)))+cos(s_f-s(ne(m, 2)))+cos(s_f-s(ne(m, 3)))+cos(s_f-s(ne(m, 4))))
            ds = new_energy - old_energy
            ! 上で計算したエネルギーよりも生成した一様乱数の値が小さければスピン反転を行う
            if(rnd2(i) < exp(-beta*dble(ds))) then
                ! 反転したスピンの数を数える
                nhit = nhit + 1
                ! スピンの反転を行う
                s(m) = s_f
            endif
        enddo
        ! 反転させたスピンの数が占める全体に対しての割合を求める
        hit = dble(nhit) / dble(n)
    end subroutine ising

    ! 角度etaを中心軸としたthetaの対称を取る
    double precision function flip(theta, eta)
        double precision,intent(in) :: theta, eta
        flip = mod(theta + 2.0d0 * (eta - theta), 2.0d0*pi)
    end function

    ! magnet    : 磁化を計算する
    ! mgn       : ひとつの格子数が平均的に持つ磁化を調べる
    ! s         : 格子点が持つスピンの向きを格納する配列
    ! n         : すべての格子点の数
    subroutine magnet(mgn, s, n)
        integer, intent(in) :: n
        double precision, intent(in) :: s(n)
        double precision, intent(out) :: mgn
        integer i
        double precision m_x, m_y, m
        
        ! 全体のスピン配位を足し合わせることによって全体の磁化を調べる
        m_x = 0.0d0
        m_y = 0.0d0

        do i=1,n
            m_x = m_x + cos(s(i))
            m_y = m_y + sin(s(i))
        enddo
        m = sqrt(m_x*m_x + m_y*m_y)
        ! 全体の磁化の割合を調べる
        mgn = m / dble(n)
    end subroutine magnet

    subroutine hamiltonian(energy, s, n, ne)
        integer, intent(in) :: n, ne(n,4)
        double precision, intent(in) :: s(n)
        double precision, intent(out) :: energy
        double precision delta
        integer i
        energy = 0.0d0
        delta = 0.0d0
        do i = 1, n
            delta = -1.0d0 * (cos(s(i) - s(ne(i, 1))) + cos(s(i) - s(ne(i, 2))) + cos(s(i) - s(ne(i, 3))) + cos(s(i) - s(ne(i, 4))))
            energy = energy + delta
        end do
        energy = energy / dble(n)
    end subroutine hamiltonian
      
    ! リストベクトル：隣接する格子点の番地を計算
    ! ne(m,direction) 
    ! m : 番地番号
    ! direction : 1->右隣, 2->上隣, 3->左隣, 4->下隣
    ! n1 : x軸方向の一周分の長さ
    ! n2 : y軸方向の一周分の長さ
    subroutine list_vector(ne, n1, n2)
        integer, intent(in) :: n1, n2
        integer, intent(out) :: ne((n1*n2),4)
        integer m1, m2, i1, i2

        ! i2=y:1 ~ n2
        do i2 = 1, n2
            ! m2:(y-1)*n1 = (y-1)L
            ! m2 is pivot of y
            m2 = (i2-1)*n1
            ! i1=x:1 ~ n1
            do i1 = 1, n1
                ! m1:(y-1)L+x
                m1 = m2 + i1
                ! ne(m, 1) = {(y-1)L + x%n1} + 1
                ! (一番左の格子点 + x)の右隣の格子点
                ! => m1の右隣(x軸正)の格子点
                ne(m1,1) = ( m2 + mod(i1,n1) ) + 1

                ! (y%n2 * n1) + x
                ! (捜索している格子点の一番右の格子点) + x
                ! => m1の上隣(y軸正)の格子点
                ne(m1,2) = mod(i2,n2)*n1 + i1

                ! ( (一番右の格子点) + L + (x - 2) )の格子点
                ! x=1の時 : L + m2 i.e. m2の一番右の格子点 (周期境界条件より正しい)
                ! x=2の時 : m2 + 1 i.e. m2 + 2の一個前の格子点
                ! x=>3の時 : m2 + x - 1 i.e. m2 + xの一個前の格子点
                ! => m1の左隣(z軸負)の格子点
                ne(m1,3) = m2 + mod(i1-2+n1,n1) + 1

                ! ( (y-2)の一番右の格子点 ) + x
                ! => m1の下隣(y軸負)の格子点
                ne(m1,4) = mod(i2-2+n2,n2)*n1 + i1
            enddo
        enddo
    end subroutine list_vector
end module subprog

program ising2d
    use extension
    use subprog
    implicit none
    double precision, allocatable :: s(:)
    integer, allocatable :: ne(:,:)
    double precision :: beta, hit, mgn, thit, tmgn, tenergy1, tenergy2, energy
    double precision :: specific_heat
    integer :: n1=32, n2=32, nbeta=180
    integer :: init=0, therm=500, inter=10, niter=5000
    integer :: n, i, j, nb
    open(10, file='output.txt')
    n=n1*n2

    allocate (s(n), ne(n,4))
    ! 毎回異なるランダム数を生成するためにSeedを設定する
    call random_seed
    ! リストベクトルを実行することによって、隣り合う格子点の番地を配列neに格納する
    call list_vector(ne, n1, n2)
    
    ! 設定条件の表示
    ! x軸方向の格子数, y軸方向の格子数, 初期条件, 物量の計算をするための回数, 磁化を計算する配位数
    write(*,200) n1, n2, init, therm, inter, niter

200 format('size:',2I3,' hot(0)/cold(1):',I2, &
    ' Therm:',I5,' Int:',I5,' NConf:',I5)
    
    ! それぞれの温度設定で計算する
    do nb=0,nbeta
        ! beta = 1 / T
        beta=0.2+nb*0.01
        ! スピンの初期化を行う
        call inits(0, s, n)
        ! thit = Total hits
        ! tmgn = Total magnetization
        thit=0.0
        tmgn=0.0
        tenergy1 = 0.0d0
        tenergy2 = 0.0d0
        ! スピンの反転をボルツマンの重みに従って行う
        ! 先にスピンの反転を行うことによってそれぞれの影響を考えた結果を得ることを期待する
        ! 最初のこの操作は隣り合う配位が正しい分布になるように行っている
        ! この最初のising()の結果は平衡状態になる前の状態なので捨てる
        do i=1,therm
            call ising(beta, s, n, hit, ne)
        enddo
        
        ! 物理量を求める操作をniter回行い平均を取ることで釣り合いを取る。
        do i=1,niter
            ! 磁化を計算する前にスピンの反転をinter回行う
            do j=1,inter
                call ising(beta, s, n, hit, ne)
                ! 反転したスピンの割合の総和を求める
                thit=thit+hit
            enddo
            ! 磁化の計算を行う
            call magnet(mgn, s, n)
            ! 磁化の合計にmgnを足し合わせる
            tmgn=tmgn+sqrt(mgn**2)
            ! エネルギーの計算を行う
            call hamiltonian(energy, s, n, ne)
            tenergy1 = tenergy1 + energy
            tenergy2 = tenergy2 + energy*energy
        enddo
        ! これまでの反転したスピンの割合の総和をスピン反転手続きを実行した数で割ることにより全体に対しての反転したスピンの割合を調べる
        tmgn=tmgn / dble(niter)
        thit=thit / dble(inter*niter)
        tenergy1 = tenergy1 / dble(niter)
        tenergy2 = tenergy2 / dble(niter)

        specific_heat = n * (tenergy2 - tenergy1*tenergy1) * beta * beta
        
        call pbflush()
        write(*,300) beta, tmgn, thit, tenergy1, specific_heat
        call pbout(nb, nbeta, .true.)
        write(10, *) 1/beta, tmgn, thit, tenergy1, specific_heat
    enddo
    close(10)
300 format('beta=',f8.5,' mgn=',f9.6,' hit=',f9.5,' energy=',f0.5,' specific heat=', f0.5)
end program ising2d