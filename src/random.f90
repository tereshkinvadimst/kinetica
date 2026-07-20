module kinetica_random
    use kinetica_kinds
    implicit none
    private

    public::poisson
    public::gauss

    
contains

    function poisson(lambda) result(k)
        !> Среднее распределение
        real(f64),intent(in),value::lambda
        !> Случайное число с плотностью вероятности Пуассона
        integer(i32)::k
        real(f64)::l
        real(f64)::p
        real(f64)::u

        l = exp(-lambda)
        k = -1
        p = 1._f64
        do
            k = k + 1
            call random_number(u)
            p = p * u
            if (p <= L) exit
        end do
    end function poisson

    !> Генерация числа в Гауссовском распределении по методу полярных координат
    function gauss() result(z)
        real(f64)::z
        real(f64)::u1, u2, s
        do
            call random_number(u1)
            call random_number(u2)
            u1 = 2._f64 * u1 - 1._f64
            u2 = 2._f64 * u2 - 1._f64
            s = u1*u1 + u2*u2
            if (s > 0._f64 .and. s < 1._f64) exit
        end do
        z = u1 * sqrt(-2._f64 * log(s) / s)
  end function gauss
end module kinetica_random