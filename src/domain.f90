module kinetica_domain
    use kinetica_kinds
    implicit none
    private

    real(f64),public::box_min_position(3)       !< Минимальное положение расчётной области
    real(f64),public::box_max_position(3)       !< Максимальное положение расчётной области
    real(f64),public::box_size(3)               !< Размер расчётной области
    integer(i32),public::boundary_conditions(6) !< Тип граничных условий
    ! Зарезервированные константы для сторон куба
    integer(i32),parameter,public::i_left  = 1  !< Левая сторона куба
    integer(i32),parameter,public::i_right = 2  !< Правая сторона куба
    integer(i32),parameter,public::i_top   = 3  !< Верхняя сторона куба
    integer(i32),parameter,public::i_down  = 4  !< Нижняя сторона куба
    integer(i32),parameter,public::i_front = 5  !< Передняя сторона куба
    integer(i32),parameter,public::i_back  = 6  !< Задняя сторона куба

end module kinetica_domain