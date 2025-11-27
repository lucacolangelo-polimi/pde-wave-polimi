#ifndef SHAPEFUNCTIONS_HPP
#define SHAPEFUNCTIONS_HPP

#include <Eigen/Dense>
#include "Mesh.hpp"

namespace ShapeFunctions {

    // Calculate the local stiffness matrix (3x3) for element e
    // K_loc = integrale(grad(phi_i) . grad(phi_j))
    Eigen::Matrix3d computeLocalStiffness(size_t e, const Mesh& mesh);

    // Calculate the local mass matrix (3x3) for element e
    // M_loc = integrale(phi_i * phi_j)
    Eigen::Matrix3d computeLocalMass(size_t e, const Mesh& mesh);

}

#endif