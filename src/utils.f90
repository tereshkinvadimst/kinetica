module kinetica_utils
    use kinetica_kinds
    implicit none
    private

    public::resize

    interface resize
        module procedure resize_matrix_f64
        module procedure resize_vector_f64
        module procedure resize_vector_i32
    end interface resize
    
contains

    subroutine resize_matrix_f64(matrix, new_capacity)
        real(f64),allocatable,intent(inout)::matrix(:,:)
        integer(i32),intent(in)::new_capacity(2)
        real(f64),allocatable::tmp(:,:)

        allocate(tmp(new_capacity(1), new_capacity(2)))
        if(allocated(matrix)) tmp(1 : size(matrix, dim=1), 1 : size(matrix,dim=2)) = matrix
        call move_alloc(from=tmp, to=matrix)
    end subroutine resize_matrix_f64

    
    subroutine resize_vector_f64(vector, new_capacity)
        real(f64),allocatable,intent(inout)::vector(:)
        integer(i32),intent(in),value::new_capacity
        real(f64),allocatable::tmp(:)

        allocate(tmp(new_capacity))
        if(allocated(vector)) tmp(1 : size(vector)) = vector
        call move_alloc(from=tmp, to=vector)
    end subroutine resize_vector_f64
    
    subroutine resize_vector_i32(vector, new_capacity)
        integer(i32),allocatable,intent(inout)::vector(:)
        integer(i32),intent(in),value::new_capacity
        integer(i32),allocatable::tmp(:)

        allocate(tmp(new_capacity))
        if(allocated(vector)) tmp(1 : size(vector)) = vector
        call move_alloc(from=tmp, to=vector)
    end subroutine resize_vector_i32

end module kinetica_utils