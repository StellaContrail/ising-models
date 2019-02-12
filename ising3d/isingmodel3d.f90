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

module ising3d
    implicit none
contains
    ! Initialize all spins
    ! s : spins array
    ! mode : 1 or -1 => Set all spins 1/-1, other => random on each spin
    subroutine init(s, mode)
        double precision,intent(out) :: s(:)
        integer,intent(in) :: mode
        double precision,allocatable :: rand(:)
        integer i, n
        n = size(s)
        allocate(rand(n))

        if (mode == -1 .or. mode == 1) then
            s(:) = dble(mode)
        else
            call random_number(rand)
            do i = 1, n
                if (rand(i) > 0.5d0) then
                    s(i) = 1.0d0
                else
                    s(i) = -1.0d0
                end if
            end do
        end if

        deallocate(rand)
    end subroutine

    subroutine list_vector(ne, n1, n2, n3)
        integer,intent(out) :: ne(:,-3:)
        integer,intent(in) :: n1, n2, n3
        integer m1, m2, m3, i1, i2, i3

        ! それぞれのmはそれぞれの軸についての行/列のpivotを表す
        do i3 = 1, n3
            m3 = (i3-1)*n1*n2
            do i2 = 1, n2
                m2 = m3 + (i2-1)*n1
                do i1 = 1, n1
                    m1 = m2 + i1

                    ! 周期境界条件を考えなければいけないのはn1 x n2 x n3の立方体の面にある格子点についてのみ
                    ! すなわち、それぞれのx座標、y座標、z座標が1またはn1/n2/n3を取るとき
                    if (i1 == n1) then
                        ne(m1, 1) = m2 + 1
                    else
                        ne(m1, 1) = m2 + i1 + 1
                    end if

                    if (i1 == 1) then
                        ne(m1, -1) = m2 + n1
                    else
                        ne(m1, -1) = m2 + i1 - 1
                    end if

                    if (i2 == n2) then
                        ne(m1, 2) = m3 + i1
                    else
                        ne(m1, 2) = m1 + n1
                    end if

                    if (i2 == 1) then
                        ne(m1, -2) = m3 + (n2-1)*n1 + i1
                    else
                        ne(m1, -2) = m1 - n1
                    end if

                    if (i3 == n3) then
                        ne(m1, 3) = (i2-1)*n1 + i1
                    else
                        ne(m1, 3) = m1 + n1*n2
                    end if

                    if (i3 == 1) then
                        ne(m1, -3) = (n3-1)*n1*n2 + (i2-1)*n1 + i1
                    else
                        ne(m1, -3) = m1 - n1*n2
                    end if
                end do
            end do
        end do
    end subroutine

    subroutine ising(beta, s, hitrate, ne)
        double precision,intent(in) :: beta
        double precision,intent(inout) :: s(:)
        double precision,intent(out) :: hitrate
        integer, intent(in) :: ne(:, -3:)
        integer i, n, m, hits
        double precision ds
        double precision,allocatable :: rand1(:), rand2(:)
        ! FORTRANではsubroutineやfunction内で定義時に初期化した変数は再び呼ばれても初期化されず、前に使われたものを再利用される。
        ! そのため別途このように初期化する必要あり。
        m = 0
        hits = 0
        n = size(s)
        allocate(rand1(n), rand2(n))

        call random_number(rand1)
        call random_number(rand2)

        do i = 1, n
            ! Choose a random point
            ! Increase m by 1 since m = 0 doesn't exist
            m = int(n*rand1(i)) + 1
            if (m > n) then
                m = n
            end if
            ds = 2*s(m)*(s(ne(m,1))+s(ne(m,-1))+s(ne(m,2))+s(ne(m,-2))+s(ne(m,-3))+s(ne(m,3)))
            if (rand2(i) < exp(-beta*dble(ds))) then
                hits = hits + 1
                s(m) = -1.0d0 * s(m)
            end if
        end do

        hitrate = dble(hits)/dble(n)
        deallocate(rand1, rand2)
    end subroutine

    subroutine magnet(mgn, s)
        double precision,intent(in) :: s(:)
        double precision,intent(out) :: mgn
        integer i, n
        double precision m
        m = 0
        n = size(s)
        do i = 1, n
            m = m + s(i)
        end do
        mgn = m / dble(n)
    end subroutine

    subroutine hamiltonian(energy, s, ne)
        double precision,intent(out) :: energy
        double precision,intent(in) :: s(:)
        integer,intent(in) :: ne(:, -3:)
        integer n, i
        double precision delta
        n = size(s)
        do i = 1, n
            delta = -s(i)*(s(ne(i,-3))+s(ne(i,-2))+s(ne(i,-1))+s(ne(i,1))+s(ne(i,2))+s(ne(i,3)))
            energy = energy + delta
        end do
        energy = energy / dble(n)
    end subroutine

    subroutine debug(s, n1, n2, n3)
        implicit none
        double precision,intent(in) :: s(:)
        integer,intent(in) :: n1, n2, n3
        integer z, y, x

        do z = 1, n3
            write (*, '(a, i0)') "Z=", z
            do y = n2, 1, -1
                do x = 1, n1
                    write (*, '(3(F4.1, x))', advance='no') s((z-1)*9+(y-1)*3+x)
                end do
                write (*, *)
            end do
            write (*, *)
        end do
        write (*, *)
    end subroutine
end module

program main
    use extension
    use ising3d
    implicit none
    double precision,allocatable :: s(:)
    integer,allocatable :: ne(:, :)
    double precision beta, hitrate, mgn, hit_total, mgn_total, energy, t_energy1, t_energy2, specific_heat
    integer :: n1 = 15, n2 = 15, n3 = 15, nbeta_max = 50, mode = 0
    integer :: trashcount = 1500, waitcount = 30, loopcount = 15000 ! the first plotted data was produced using all values divided by three
    integer :: n, i, j, nbeta
    integer,parameter :: fo = 10
    open(10, file='output.txt')
    n = n1*n2*n3
    allocate(s(n), ne(n, -3:3))
    call list_vector(ne, n1, n2, n3)

    write (*, 100) n1, n2, n3, mode, trashcount, waitcount, loopcount

    do nbeta = 0, nbeta_max
        beta = 0.1d0 + nbeta*0.01d0
        call init(s, mode)

        hit_total = 0.0d0
        mgn_total = 0.0d0

        do i = 1, trashcount
            call ising(beta, s, hitrate, ne)
        end do
        ! 物理量の釣り合いをするためのLoop
        do i = 1, loopcount
            do j = 1, waitcount ! Ising Fieldの釣り合いをするためのLoop
                call ising(beta, s, hitrate, ne)
                hit_total = hit_total + hitrate
            end do
            call magnet(mgn, s)
            mgn_total = mgn_total + abs(mgn)
            call hamiltonian(energy, s, ne)
            t_energy1 = t_energy1 + energy
            t_energy2 = t_energy2 + energy * energy
        end do
        mgn_total = mgn_total / dble(loopcount)
        hit_total = hit_total / dble(loopcount*waitcount)
        t_energy1 = t_energy1 / dble(loopcount)
        t_energy2 = t_energy2 / dble(loopcount)

        specific_heat = n * (t_energy2 - t_energy1*t_energy1) * beta * beta
    
        !call debug(s, n1, n2, n3)
        call pbflush()
        write (*, 200) beta, mgn_total, hit_total, t_energy1, specific_heat
        call pbout(nbeta, nbeta_max, .true.)
        write (fo, *) beta, 1/beta, mgn_total, hit_total, t_energy1, specific_heat
    end do
    close(fo)
    deallocate(s, ne)
100 format('size:', I0, 'x', I0, 'x', I0, ' mode=', I0, ' TrashCount:', I0, ' WaitCount:', I0, ' LoopCount:', I0)
200 format('beta=', F8.5, ' Magnetization:', F9.6, ' HitRate:', F0.5, ' Energy:', F0.5, ' SpecificHeat:', F0.5)
end program