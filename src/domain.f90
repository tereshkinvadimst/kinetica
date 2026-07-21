module kinetica_domain
    use kinetica_kinds
    implicit none
    private

    real(f64),public::box_min_position(3)       !< Минимальное положение расчётной области
    real(f64),public::box_max_position(3)       !< Максимальное положение расчётной области
    real(f64),public::box_size(3)               !< Размер расчётной области
    logical,public::periodic(3)                 !< Периодические граничные условия по осям
    logical,public::is_wall(6)                  !< Является ли граница стенкой
    ! Зарезервированные константы для сторон куба
    integer(i32),parameter,public::i_left  = 1  !< Левая сторона куба
    integer(i32),parameter,public::i_right = 2  !< Правая сторона куба
    integer(i32),parameter,public::i_top   = 3  !< Верхняя сторона куба
    integer(i32),parameter,public::i_bot   = 4  !< Нижняя сторона куба
    integer(i32),parameter,public::i_front = 5  !< Передняя сторона куба
    integer(i32),parameter,public::i_back  = 6  !< Задняя сторона куба
    ! Зарезервированные константы для сопостовления осей
    integer(i32),parameter,public::axis_min_face(3) = [i_left,  i_bot, i_back ]
    integer(i32),parameter,public::axis_max_face(3) = [i_right, i_top, i_front]

end module kinetica_domain