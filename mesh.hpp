#ifndef MESH_HPP
#define MESH_HPP

/*
The file defines the fundamental data structures for representing a 2D triangular mesh for the finite element method.
*/

#include <vector>
#include <array>
#include <string>    
#include <cstddef>   

struct Node {
    int id;                     // unique identifier for the node
    double x, y;
};

struct Element {
    int id;                      // unique identifier for the element
    std::array<int,3> nodes;     // 3 vertices indices, every element is a triangle
};

struct ElementGeometryData {
    double area;                        // area of the triangle
    std::array<double, 2> grad_phi1;    // gradient of the basis function associated to vertex 1
    std::array<double, 2> grad_phi2;
    std::array<double, 2> grad_phi3;
};

class Mesh {
public:
    std::vector<Node> nodes;        // list of all nodes in the mesh, Each node has a spatial coordinate(x,y)
    std::vector<Element> elements;  // list of all elements in the mesh
    
    // List of global indices of boundary nodes (useful for applying Dirichlet BCs)
    std::vector<int> boundary_node_indices; 

    // List to store pre-calculated geometric data, once per element. (in this way we avoid redundant calculations during assembly)
    std::vector<ElementGeometryData> geometry_data; 
    
    // Getters
    size_t get_num_elements() const { 
        return elements.size(); 
    }
    
    // Get number of nodes
    size_t get_num_nodes() const {
        return nodes.size();
    }

    // Methods
    void generate_simple_mesh();            // Generate a simple predefined mesh (for testing)
    void read_from_file(const std::string& filename);  // Read mesh from a file (e.g., in a simple custom format)
    void find_boundary_nodes(double tolerance = 1e-9);  // Identify boundary nodes based on coordinates
    void compute_all_element_geometry();    // Pre-compute geometric data for all elements. 
};

#endif
