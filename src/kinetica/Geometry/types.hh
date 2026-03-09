#ifndef MF_KINETICA_GEOMETRY_TYPES_H
#define MF_KINETICA_GEOMETRY_TYPES_H
#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace mf {

using Vector3 = Eigen::Vector3d;
using PlaneT  = Eigen::Hyperplane<double, 3>;

}  // namespace mf

#endif  // MF_KINETICA_GEOMETRY_TYPES_H