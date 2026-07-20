module kinetica_constants
    use kinetica_kinds
    implicit none
    private
    !> Постоянная Авогадро
    real(f64),parameter,public::n_a    = 6.02214076e23_f64
    !> Постоянная Больцмана
    real(f64),parameter,public::k_b    = 1.380649e-23_f64
    !> Универсальная газовая постоянная
    real(f64),parameter,public::r_gas  = k_b * n_a
    !> Число пи
    real(f64),parameter,public::pi     = 3.14159265358979323846_f64
    !> Обратное число пи = 1 / пи
    real(f64),parameter,public::inv_pi = 1 / pi
end module kinetica_constants