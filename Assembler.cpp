#include "Assembler.hpp"
#include <iostream>

Assembler::Assembler(const Mesh& mesh_input) : mesh(mesh_input) {}

void Assembler::Assemble(SparseMatrix& K, Vector& M) {
    
    size_t num_nodes = mesh.get_num_nodes();    //setup
    size_t num_elements = mesh.get_num_elements();

    std::cout << "Starting Assembly (Nodes: " << num_nodes 
              << ", Elements: " << num_elements << ")..." << std::endl;

    // mass is initialized to zero
    M.resize(num_nodes);
    M.setZero();

    // Triplets list
    std::vector<Triplet> triplets;

    //We reserve memory: each triangular element has 3 nodes, so 3x3=9 local interactions
    triplets.reserve(num_elements * 9); 

    for (size_t e = 0; e < num_elements; ++e) {
        
        const auto& elem = mesh.elements[e];
        const auto& geo = mesh.geometry_data[e];

        double area = geo.area;

        const std::array<double, 2>* grads[3] = {
            &geo.grad_phi1, 
            &geo.grad_phi2, 
            &geo.grad_phi3
        };

        //mass lumped contribution
        double lumped_mass_contribution = area / 3.0;

        for (int i = 0; i < 3; ++i) {
            int global_node_id = elem.nodes[i];
            M[global_node_id] += lumped_mass_contribution;
        }

        // K_ij = Area * (∇phi_i · ∇phi_j)
        
        for (int i = 0; i < 3; ++i) {       
            for (int j = 0; j < 3; ++j) {    
                
                double grad_dot_prod = (*grads[i])[0] * (*grads[j])[0] + 
                                       (*grads[i])[1] * (*grads[j])[1];
                
                double value = area * grad_dot_prod;

                int global_row = elem.nodes[i];
                int global_col = elem.nodes[j];

                triplets.push_back(Triplet(global_row, global_col, value));
            }
        }
    }

    K.resize(num_nodes, num_nodes);
    K.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << "Assembly Completed." << std::endl;
}