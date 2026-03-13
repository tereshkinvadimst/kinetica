#include "kinetica/DSMC/collider.hh"

#include <random>

auto mf::collider(Particles&                            particles,
                  std::span<const Particles::size_type> particles_id,
                  Particles::value_type                 volume,
                  Particles::value_type                 dt,
                  Particles::value_type                 sigma_g_max,
                  random&                               gen) -> Particles::value_type {
    using value_type = Particles::value_type;
    using size_type  = Particles::size_type;

    static std::uniform_real_distribution<value_type> u01(value_type{}, value_type{1});
    static std::normal_distribution<value_type>       n01(value_type{}, value_type{1});

    // Новое сечение столкновения
    value_type new_sigma_g = sigma_g_max;

    // Число частиц, которые необходимо между собой "столкнуть"
    const auto Np = particles_id.size();
    // Создаем равномерное распределение индексов частиц
    std::uniform_int_distribution<size_type> rid(0, Np - 1);
    // Если система слишком маленькая, то пропускаем её
    if (Np < 2) return new_sigma_g;
    // Вычисляем интенсивность столкновений
    const auto collision_rate = sigma_g_max / volume;
    // Примитивная оценка среднего числа частиц в ячейке
    const auto Navg = Np;
    // Вычисляем среднее число столкновений в ячейке
    const auto Ncoll_avg = 0.5 * Np * Navg * particles.W * collision_rate * dt;
    // Создаём распределение Пуассона
    std::poisson_distribution<size_type> poisson(Ncoll_avg);
    // Вычисляем число сталкивающихся частиц из распределения Пуассона
    const auto Ncoll = poisson(gen);
    for (size_type collision{}; collision < Ncoll; ++collision) {  // Чикл по числу коллизий
        //-------------------------------------------------------------------------------------------
        // Получаем случайный индекс первой частицы
        const auto i = particles_id[rid(gen)];
        // Получаем случайный индекс второй частиц (и проверяем, чтобы он не равнялся индексу первой)
        auto tmp = particles_id[rid(gen)];
        while (tmp == i) {
            tmp = particles_id[rid(gen)];
        }
        const auto j = tmp;
        //-------------------------------------------------------------------------------------------
        // Вычисляем разность скоростей между частицами
        const auto dux = particles.ux[i] - particles.ux[j];
        const auto duy = particles.uy[i] - particles.uy[j];
        const auto duz = particles.uz[i] - particles.uz[j];
        // Вычисляем квадрат относительной скорости
        const auto gij2 = dux * dux + duy * duy + duz * duz;
        // Вычисляем относительную скорость
        const auto gij = std::sqrt(gij2);
        // Вычисляем половину относительной скорости
        const auto hgij = value_type{0.5} * gij;
        // Вычисляем сечение столкновения
        const auto sigma_g = particles.sigma(gij) * gij;
        // Проверяем, что вычисленное сечение столкновения больше максимального
        if (sigma_g > sigma_g_max) new_sigma_g = sigma_g;
        if (u01(gen) < sigma_g / new_sigma_g) {  // Принимаем столкновения с равновероятной вероятностью
            //-------------------------------------------------------------------------
            // Считаем скорость центра масс сталкивающихся частиц
            const auto ucm_x = 0.5 * (particles.ux[i] + particles.ux[j]);
            const auto ucm_y = 0.5 * (particles.uy[i] + particles.uy[j]);
            const auto ucm_z = 0.5 * (particles.uz[i] + particles.uz[j]);
            //-------------------------------------------------------------------------
            // Выбираем случайного направление после столкновения
            auto nx = n01(gen);
            auto ny = n01(gen);
            auto nz = n01(gen);
            // Вычисляем норму случайного направления
            const auto norm = std::sqrt(nx * nx + ny * ny + nz * nz);
            // Нормируем направление
            nx /= norm;
            ny /= norm;
            nz /= norm;
            //-------------------------------------------------------------------------
            /// Обновляем скорости
            particles.ux[i] = ucm_x + hgij * nx;
            particles.uy[i] = ucm_y + hgij * ny;
            particles.uz[i] = ucm_z + hgij * nz;

            particles.ux[j] = ucm_x - hgij * nx;
            particles.uy[j] = ucm_y - hgij * ny;
            particles.uz[j] = ucm_z - hgij * nz;
        }
    }
    return new_sigma_g;
}