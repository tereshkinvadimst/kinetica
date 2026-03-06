#ifndef MF_KINETICA_STATS_H
#define MF_KINETICA_STATS_H
#pragma once
#include <string>
#include <unordered_map>

namespace mf {

class Stats final {
   public:
    using value_type = double;
    using size_type  = std::size_t;

    Stats()          = default;

    void initStat(std::string key);
    void addStat(std::string key, value_type value);
    void printHeader() const;
    void printStats(value_type time) const;

   private:
    std::unordered_map<std::string, value_type> data_;
    std::unordered_map<std::string, size_type>  n_writes_;
};

}  // namespace mf

#endif  // MF_KINETICA_STATS_H