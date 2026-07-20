program kinetica
    use kinetica_kinds
    use kinetica_random
    use kinetica_particles
    use kinetica_species
    use kinetica_xyz
    implicit none

    call random_init(repeatable=.false., image_distinct=.true.)

    call add_species("Ar", 0.03995_f64, 0.81_f64, 1.4_f64, 4.17e-10_f64, 273._f64)

    call generate_particles(1._f64, 1e26_f64, 300._f64, [0._f64, 0._f64, 0._f64], 1, [0._f64, 0._f64, 0._f64], [100e-9_f64, 100e-9_f64, 100e-9_f64])

    call write_xyz("test")
end program kinetica
