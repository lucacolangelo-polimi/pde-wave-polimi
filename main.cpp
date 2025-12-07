#include "Mesh.hpp"
#include "ShapeFunctions.hpp"
#include <iostream>
#include <cmath>

int main() {
    Mesh my_mesh;
    
    // generating a simple test mesh
    my_mesh.generate_simple_mesh(); 
    my_mesh.find_boundary_nodes();

    // Display boundary nodes
    my_mesh.compute_all_element_geometry();


if (my_mesh.get_num_elements() > 0) {
        size_t e = 0; 
        const auto& geo = my_mesh.geometry_data[e]; 
        
        std::cout << "\n--- TEST GEOMETRIA ELEMENTO " << e << " ---" << std::endl;
        std::cout << "Area: " << geo.area << std::endl;
        std::cout << "Grad(phi1): (" << geo.grad_phi1[0] << ", " << geo.grad_phi1[1] << ")" << std::endl;
        std::cout << "Grad(phi2): (" << geo.grad_phi2[0] << ", " << geo.grad_phi2[1] << ")" << std::endl;

        // Eigen::Matrix3d it's correct
        Eigen::Matrix3d K = ShapeFunctions::computeLocalStiffness(e, my_mesh);
        
        std::cout << "\nStifness matrix (K):" << std::endl;
        std::cout << K << std::endl; 

        // Mass Lumped
        Eigen::Matrix3d M = ShapeFunctions::computeLocalMass(e, my_mesh);
        
        std::cout << "\n Mass Matrix Lumped (M):" << std::endl;
        std::cout << M << std::endl;
        
        std::cout << "\n------" << std::endl;
        
        //Mass matrix should be diagonal
        bool isDiagonal = (std::abs(M(0,1)) < 1e-9 && std::abs(M(0,2)) < 1e-9);
        std::cout << "Diagonal matrix? " << (isDiagonal ? "SI (Ok)" : "NO (Errore)") << std::endl;

        // Stiffness matrix row sums should be ~0
        double rowSum0 = K(0,0) + K(0,1) + K(0,2);
        std::cout << "rows sum  K (needs to be ~0): " << rowSum0 << std::endl;
    }


    return 0;
}