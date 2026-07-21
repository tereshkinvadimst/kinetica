module kinetica_cells
    use kinetica_kinds
    implicit none
    private

    integer(i32),protected,public::n_cells(3) !< Число ячеек по каждой оси
    real(f64),protected,public::cell_sigmag_max(:,:,:)
    
contains
    
end module kinetica_cells