
#include "lab/Properties/stats.hh"

#include <iomanip>
#include <iostream>
#include <ranges>

namespace ranges          = std::ranges;
namespace views           = ranges::views;

constexpr auto streamsize = 20;

void           mf::Stats::initStat(std::string key) {
    data_[key]     = value_type{};
    n_writes_[key] = size_type{};
}
void mf::Stats::addStat(std::string key, value_type value) {
    data_[key] += value;
    ++n_writes_[key];
}
void mf::Stats::printHeader() const {
    std::cout << std::left << std::setw(streamsize) << "time [s]";
    for (const auto& key : data_ | views::keys) {
        std::cout << std::right << std::setw(streamsize) << key;
    }
    std::cout << '\n';
}
void mf::Stats::printStats(value_type time) const {
    constexpr auto pressicion = 8;
    std::cout << std::left << std::setw(streamsize) << std::scientific << time;
    for (const auto& [key, value] : data_) {
        std::cout << std::right << std::setw(streamsize) << std::setprecision(pressicion) << value / n_writes_.at(key);
    }
    std::cout << '\n';
}