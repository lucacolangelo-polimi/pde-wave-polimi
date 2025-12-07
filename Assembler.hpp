#ifndef ASSEMBLER_HPP
#define ASSEMBLER_HPP

#include "Mesh.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>

class Assembler {
public:
 
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using Vector       = Eigen::VectorXd;
    using Triplet      = Eigen::Triplet<double>; 

    Assembler(const Mesh& mesh_input);
    //main function to assemble global stiffness matrix K and mass vector M
    void Assemble(SparseMatrix& K, Vector& M);

private:
    const Mesh& mesh; 
};

#endif 