#include "Mesh.h"
#include <fstream>
#include <iostream>
#include <iomanip> //set precision




// XYZ //
void Mesh::load_vertices_xyz(std::string filename) {
    // X Y Z file
    double x, y, z;
    unsigned int vertice_id = 0;
    std::ifstream src(filename);
    while (src >> x >> y >> z) {
        vertice_id += 1;
        Mesh::add_vertice(x,y,z,vertice_id);
    }
    src.close();
}
void Mesh::save_vertices_xyz(std::string filename) {
    std::ofstream src(filename);
    src << std::scientific;
    src << std::setprecision(8);
    for (Vertice* v : vertices) {
        src << v->coor->x << ' ';
        src << v->coor->y << ' ';
        src << v->coor->z << '\n';
    }
    src.close();
}
void Mesh::load_elements_with_vertice_ids(std::string filename) {
    //each line one element, and list of id separated by a space delimiter
    std::ifstream src(filename);
    std::string delimiter = " ";
    size_t pos = 0;
    std::string token;
    for( std::string line; getline( src, line ); ) {
        //one line one element
        std::vector<unsigned int> ids;
        while ((pos = line.find(delimiter)) != std::string::npos) {
            token = line.substr(0, pos);
            line.erase(0, pos + delimiter.length());
            ids.push_back(std::stoi(token));
        }
        add_element(ids);
    }
    src.close();
}


// TetGen //
void Mesh::load_vertices_tetgen(std::string filename) {
    // X Y Z file
    double x, y, z;
    unsigned int vertice_id;
    std::ifstream src(filename);
    while (src >> vertice_id >> x >> y >> z) {
        add_vertice(x,y,z,vertice_id);
    }
    src.close();
}
void Mesh::save_vertices_tetgen(std::string filename) {
    std::ofstream src(filename);
    unsigned int count = 1;
    for (Vertice* v : vertices) {
        src << std::fixed;
        src << count << ' ';
        src << std::scientific;
        src << std::setprecision(8);
        src << v->coor->x << ' ';
        src << v->coor->y << ' ';
        src << v->coor->z << '\n';
        count++;
    }
    src.close();
}
void Mesh::load_elements_tetgen(std::string filename) {
    //each line one element, and list of id separated by a space delimiter
    std::ifstream src(filename);
    unsigned int read_value;
    std::vector<unsigned int> ids;
    while (src >> read_value) {
        ids.push_back(read_value);
        if (ids.size() == 5) {
            ids.erase(ids.begin());
            add_element(ids);
            ids.clear();
        }
    }
    src.close();
}


// Load from array
void Mesh::load_vertices_array(double vertices[], size_t n_vertices) {
    //load vertices from a array
    //designed for Python external call
    //vertices array must be 1D of size n_vertices*3
    for (unsigned int id=1; id!=n_vertices+1; id++) {
        add_vertice(vertices[3*id], vertices[3*id+1],
                          vertices[3*id+2], id);
    }
}

void Mesh::load_elements_array(unsigned int elements[], unsigned int type[],
                               size_t n_elements) {
    //load elements from a array
    //designed for Python external call
    //elements array must be 1D and contained the ids of elements contiguously
    //type array must be 1D with the number of vertice of the considered element
    size_t count=0;
    std::vector<unsigned int> ids;
    for (unsigned int n_to_read=0; n_to_read!=n_elements; n_to_read++) {
        ids.resize(0); //clear for reading a new elements
        for (unsigned int x=0; x!=type[n_to_read]; ) { //read ids
            ids.push_back(elements[count]);
            count++;
        }
        add_element(ids); //create the element
    }
}

