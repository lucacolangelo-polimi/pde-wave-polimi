#include "ShapeFunctions.hpp"

/**
 * @brief Computes the local mass matrix for the element.
 * * NOTE ON MASS MATRIX TYPES:
 * * 1. Consistent Mass Matrix (Standard):
 * Physically, this assumes the mass is continuously distributed across the 
 * entire triangle area. While more accurate, it results in a global matrix M 
 * that is sparse but NOT diagonal (it has off-diagonal terms). 
 * - Drawback: In explicit time-stepping (M*a = F), finding the acceleration 'a' 
 * requires solving a linear system at every time step (a = M^-1 * F). 
 * - Consequence: This is computationally expensive and slow for large meshes.
 * * 2. Lumped Mass Matrix (Diagonal):
 * This approximation "lumps" the total mass of the element onto its vertices (nodes).
 * This results in a strictly diagonal global matrix.
 * - Advantage: Inverting a diagonal matrix is trivial. We don't need to solve a 
 * linear system; acceleration is computed via simple element-wise division 
 * (a_i = F_i / M_ii).
 * - Consequence: This makes explicit time integration virtually instantaneous 
 * and extremely efficient.
 */


namespace ShapeFunctions {

    // Mass matrix
    /*Eigen::Matrix3d computeLocalMass(size_t e, const Mesh& mesh) {
        // 1. Retrieve the area pre-computed in the Mesh module
        double area = mesh.geometry_data[e].area;

        // 2. Standard formula for linear P1 triangles
        // (Area / 12) * [2 1 1; 1 2 1; 1 1 2]
        Eigen::Matrix3d M_loc;
        double coeff = area / 12.0;

        M_loc << 2, 1, 1,
                 1, 2, 1,
                 1, 1, 2;

        return M_loc * coeff;
    }
    */
    
    //LUMPED MASS MATRIX (Diagonal)
    // Useful for explicit time integration schemes
    Eigen::Matrix3d computeLocalMass(size_t e, const Mesh& mesh) {
        double area = mesh.geometry_data[e].area;
        
        // Distribute the total mass (area) equally to the 3 nodes
        double val = area / 3.0;

        return val * Eigen::Matrix3d::Identity();
    }
    

    // Stiffness matrix
    Eigen::Matrix3d computeLocalStiffness(size_t e, const Mesh& mesh) {
        // 1. Retrieve geometric data for element 'e'
        const auto& geo = mesh.geometry_data[e];
        double area = geo.area;

        // Retrieve pre-computed gradients (2-component vectors)
        // Accessing the arrays stored in the Mesh module
        const double* grads[3] = { 
            geo.grad_phi1.data(), 
            geo.grad_phi2.data(), 
            geo.grad_phi3.data() 
        };

        Eigen::Matrix3d K_loc;

        // 2. Double loop to fill the 3x3 matrix
        // K_ij = Area * (grad_phi_i * grad_phi_j)
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                
                // Dot product: (dx_i * dx_j) + (dy_i * dy_j)
                double dot_product = grads[i][0] * grads[j][0] + 
                                     grads[i][1] * grads[j][1];
                
                K_loc(i, j) = area * dot_product;
            }
        }

        return K_loc;
    }
}