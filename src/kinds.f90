module kinetica_kinds
    use iso_fortran_env, only: f32 => real32 &
                             , f64 => real64 &
                             , i32 => int32  &
                             , i64 => int64
    implicit none
    private

    public::f32
    public::f64
    public::i32
    public::i64
    
end module kinetica_kinds