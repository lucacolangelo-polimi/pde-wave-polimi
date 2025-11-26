#include "Mesh.hpp"
#include <iostream>
#include <fstream> 
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>

// A small tolerance for floating-point comparisons (needed for find_boundary_nodes)
const double TOL = 1e-9;

// Generate a simple predefined mesh (for testing)
void Mesh::generate_simple_mesh() {
    std::cout << "Test mesh generation [0, 1] x [0, 1] (2 elements P1)..." << std::endl;

    // 1. Definition Nodes (4 Nodes)
    nodes = {
        {0, 0.0, 0.0}, // Node 0 (0,0)
        {1, 1.0, 0.0}, // Node 1 (1,0)
        {2, 1.0, 1.0}, // Node 2 (1,1)
        {3, 0.0, 1.0}  // Node 3 (0,1)
    };

    // 2. Definition of Elements (2 Triangles P1)
    elements = {
        {0, {0, 1, 2}}, // Triangle 0: Nodes 0, 1, 2
        {1, {0, 2, 3}}  // Triangle 1: Nodes 0, 2, 3
    };

    std::cout << "Mesh generated with " << nodes.size() << " nodes and " 
              << elements.size() << " elements." << std::endl;
}


void Mesh::read_from_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open MSH file: " << filename << std::endl;
        return;
    }

    std::cout << "Reading mesh from MSH file: " << filename << "..." << std::endl;
    
    // Clear containers before starting
    nodes.clear();
    elements.clear();

    std::string line;
    bool reading_nodes = false;
    bool reading_elements = false;
    size_t num_entities = 0; // Total number of nodes or elements in the current section

    while (std::getline(file, line)) {
        // --- Tag Management ---
        if (line == "$Nodes") {
            reading_nodes = true;
            reading_elements = false;
            continue;
        }
        if (line == "$EndNodes") {
            reading_nodes = false;
            // Reset num_entities for the next section
            num_entities = 0;
            continue;
        }
        if (line == "$Elements") {
            reading_elements = true;
            reading_nodes = false;
            continue;
        }
        if (line == "$EndElements") {
            reading_elements = false;
            // Reset num_entities for the next section
            num_entities = 0;
            continue;
        }
        
        // Ignore headers (like $MeshFormat) or comments
        if (line.empty() || line[0] == '$' || line[0] == '#') continue;

        std::stringstream ss(line);
        
        // --- 1. Node Reading Section ---
        if (reading_nodes) {
            // The first line after $Nodes contains the total number of nodes.
            if (num_entities == 0) {
                ss >> num_entities;
                continue;
            }
            
            int id;
            double x, y, z; // GMSH always uses 3 coordinates (X, Y, Z)

            if (ss >> id >> x >> y >> z) {
                // We assume the ID corresponds to the index (0-based) for consistency.
                // We ignore the Z-coordinate for this 2D problem.
                nodes.push_back({id, x, y});
            }
            // Optimization: If you use the number of entities for reading, you could stop here.
            continue; 
        }
        
        // --- 2. Element Reading Section ---
        if (reading_elements) {
            // The first line after $Elements contains the total number of elements.
            if (num_entities == 0) {
                ss >> num_entities;
                continue;
            }
            
            // For the Msh2 format line:
            // Element_ID | Type | N. Tags | Tag 1 | Tag 2 | ... | Node ID 1 | Node ID 2 | Node ID 3
            int id, type, num_tags, n1, n2, n3;

            // Try to read the fixed part of the element line
            if (ss >> id >> type >> num_tags) {
                
                // Ignore all tags (num_tags)
                for (int i = 0; i < num_tags; ++i) {
                    int tag; 
                    ss >> tag;
                }
                
                // Element type '2' in Msh2 is a 3-node P1 triangle
                if (type == 2) { 
                    if (ss >> n1 >> n2 >> n3) {
                        // CRITICAL NOTE: GMSH IDs are 1-based. Your C++ indices are 0-based.
                        // You must convert GMSH IDs by subtracting 1 for correct array access:
                        elements.push_back({id, {n1 - 1, n2 - 1, n3 - 1}}); 
                        // The conversion n-1 is applied here!
                    }
                }
                // We ignore other types (boundary lines, tetrahedrons, etc.)
            }
            continue;
        }
    }

    std::cout << "  - Read " << nodes.size() << " nodes and " 
              << elements.size() << " elements (triangles)." << std::endl;
}

void Mesh::find_boundary_nodes(double tolerance) {
    boundary_node_indices.clear();      //The vector is emptied at the beginning of each call.
    
    // For the square [0, 1] x [0, 1], the boundary nodes are those where
    // x is close to 0 or 1, or y is close to 0 or 1.
    for (const auto& node : nodes) {
        // Check if the node is close to one of the edges
        bool on_x_boundary = (std::abs(node.x - 0.0) < tolerance) || (std::abs(node.x - 1.0) < tolerance);
        bool on_y_boundary = (std::abs(node.y - 0.0) < tolerance) || (std::abs(node.y - 1.0) < tolerance);
        
        if (on_x_boundary || on_y_boundary) {
            boundary_node_indices.push_back(node.id);
        }
    }
    
    // We remove any duplicates (even though with the sequential ID there shouldn't be any)
    std::sort(boundary_node_indices.begin(), boundary_node_indices.end());
    boundary_node_indices.erase(std::unique(boundary_node_indices.begin(), boundary_node_indices.end()), boundary_node_indices.end());

    std::cout << "Finded " << boundary_node_indices.size() 
              << " boundary nodes (Dirichlet)." << std::endl;
}


void Mesh::compute_all_element_geometry() {
    geometry_data.resize(elements.size());
    std::cout << "calculation of geometric data for " << elements.size() << " elements..." << std::endl;

    for (size_t e = 0; e < elements.size(); ++e) {
     // Obtaining the coordinates of the 3 nodes (N1, N2, N3)
     // We use the local names (n1, n2, n3) for simplicity of the formula.
        const auto& n1 = nodes[elements[e].nodes[0]];
        const auto& n2 = nodes[elements[e].nodes[1]];
        const auto& n3 = nodes[elements[e].nodes[2]];

     // 1. Components of the Jacobian (based on node 1) 
     // Gradients are computed with respect to the reference coordinate system.
     // dx_21 = x2 - x1, dy_21 = y2 - y1, etc.
        double dx_21 = n2.x - n1.x;
        double dy_21 = n2.y - n1.y;
        double dx_31 = n3.x - n1.x;
        double dy_31 = n3.y - n1.y;

        // 2. Area Calculation 
        double detJ = dx_21 * dy_31 - dx_31 * dy_21;
        
        // Check for degenerate element (area = 0)
        if (std::abs(detJ) < TOL) {
            std::cerr << "ERROR " << elements[e].id << " (Area = 0)." << std::endl;
            continue; 
        }

        double area = 0.5 * std::abs(detJ);
        geometry_data[e].area = area;

        // 3. Gradient Calculation of Basis Functions
       
        // Node 1: (d/dx, d/dy) of phi_1
        geometry_data[e].grad_phi1[0] = (dy_21 - dy_31) / detJ; // d/dx
        geometry_data[e].grad_phi1[1] = (dx_31 - dx_21) / detJ; // d/dy

        // Node 2: (d/dx, d/dy) of phi_2
        geometry_data[e].grad_phi2[0] = dy_31 / detJ; // d/dx
        geometry_data[e].grad_phi2[1] = -dx_31 / detJ; // d/dy

        // Node 3: (d/dx, d/dy) of phi_3
        geometry_data[e].grad_phi3[0] = -dy_21 / detJ; // d/dx
        geometry_data[e].grad_phi3[1] = dx_21 / detJ; // d/dy
    }

}