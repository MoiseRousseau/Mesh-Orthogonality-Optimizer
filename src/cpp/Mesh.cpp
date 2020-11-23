#include "Mesh.h"
#include <fstream>
#include <iostream>




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
        src << count << ' ';
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



