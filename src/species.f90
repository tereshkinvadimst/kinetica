module kinetica_species
    use kinetica_kinds
    implicit none
    private
    
    public::add_species

    integer(i32),protected,public::n_species = 0
    !> Имя частицы (NSpecies)
    character(len=:),allocatable,protected,public::species_names(:)
    !> Молярная масса молекулы (NSpecies)
    real(f64),allocatable,protected,public::molar_mass(:)
    !> Показатель вязкости (NSpecies)
    real(f64),allocatable,protected,public::omega(:)
    !> Параметр рассеяния (NSpecies)
    real(f64),allocatable,protected,public::alpha(:)
    !> Опорный диаметр молекулы (NSpecies)
    real(f64),allocatable,protected,public::molecule_size(:)
    !> Опорная температура молекулы (NSpecies)
    real(f64),allocatable,protected,public::reference_temperature(:)
    

contains

    subroutine add_species(name, m, omega_, alpha_, d_ref, t_ref)
        character(len=*),intent(in)::name
        real(f64),intent(in),value::m
        real(f64),intent(in),value::omega_
        real(f64),intent(in),value::alpha_
        real(f64),intent(in),value::d_ref
        real(f64),intent(in),value::t_ref


        if (n_species == 0) then
            species_names         = [name]
            molar_mass            = [m]
            omega                 = [omega_]
            alpha                 = [alpha_]
            molecule_size         = [d_ref]
            reference_temperature = [t_ref]
        else
            species_names         = [species_names, name]
            molar_mass            = [molar_mass, m]
            omega                 = [omega, omega_]
            alpha                 = [alpha, alpha_]
            molecule_size         = [molecule_size, d_ref]
            reference_temperature = [reference_temperature, t_ref]
        end if

        n_species = n_species + 1

    end subroutine add_species
    
end module kinetica_species