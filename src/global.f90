module kinetica_global
    use kinetica_kinds
    implicit none
    private

    real(f64),public::dt_coeff     !< Коэффициент уменьшения шага по времени 0.1...0.5
    real(f64),public::dl_coeff     !< Коэффициент уменьшения размера ячейки 0.3...1
    real(f64),public::dt           !< Временной шаг расчёта
    real(f64),public::max_velocity !< Максимальная скорость частицы в расчёте
    real(f64),public::lambda_min   !< Минимальная длина свободного пробега
    real(f64),public::sigmag_max   !< Максимальная интенсивность столкновения
    real(f64),public::dl           !< Размер ячейки
    real(f64),public::weight       !< Вес частиц в расчёте
    
end module kinetica_global