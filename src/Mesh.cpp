#include "Mesh.h"
#include "Connection.h"
#include "Element.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <array>
#include <iomanip> //set precision

void Mesh::decompose() {
    std::cout << "Finding mesh internal faces (decompose)" << std::endl;
    std::map<std::array<unsigned int, 4>, Connection*> unique_id_map;
    std::vector<Connection*> temp_connection;
    
    //link element to vertice instances
    for (Element* elem : elements) {
        for (unsigned int id : elem->vertice_ids) {
            elem->vertices.push_back(vertices[id-1]);
        }
    }
    for (Element* elem : elements) {
        if (elem->type == -3) { //triangle
            build_connection(0,1,-1,-1, 2, elem, unique_id_map,
                             temp_connection);
            build_connection(1,2,-1,-1, 0, elem, unique_id_map,
                             temp_connection);
            build_connection(2,0,-1,-1, 1, elem, unique_id_map,
                             temp_connection);
        }
        else if (elem->type == -4) { //quadrilateral
        
        }
        else if (elem->type == 4) { //tetrahedron
            //face 0 1 2, opposite 3
            build_connection(0,1,2,-1, 3, elem, unique_id_map,
                             temp_connection);
            //face 1 2 3, opposite 0
            build_connection(1,2,3,-1,0, elem, unique_id_map,
                             temp_connection);
            //face 0 1 3, opposite 2
            build_connection(0,1,3,-1, 2, elem, unique_id_map,
                             temp_connection);
            //face 0 2 3, opposite 1
            build_connection(0,2,3,-1,1, elem, unique_id_map,
                             temp_connection);
        }
        else if (elem->type == 5) { //pyramid
            build_connection(0,1,2,3,4, elem, unique_id_map,
                             temp_connection);
            //other are treated by tet above
        }
        else if (elem->type == 6) { //prisms
            //not optimized, so do not figure
            //TODO: I must fix the point in this case, else, they would move and degenerate the prisms. Verify
        }
        else if (elem->type == 8) { //hex
            //not optimized, so do not figure
        }
        else {
            std::cerr << std::endl;
            std::cerr << "Mesh decompose error" << std::endl;
            std::cerr << "element type not recognized: " << elem->type;
            std::cerr << std::endl;
            std::cerr << std::endl;
            exit(1);
        }
    }
    for (Connection* con : temp_connection) {
#ifdef DEBUG_MODE
        std::cout << *con << std::endl;
#endif
        if (con->vertice_up == nullptr) {
            for (Vertice* v : con->vertices) {
                v->fix_vertice();
            }
            boundary_connections.push_back(con);
        }
        else {
            connections_internal.push_back(con);
        }
    }
    n_connections_internal = connections_internal.size();
    n_connections_bc = boundary_connections.size();
    
//    for (Connection* con : connections_internal) {
//        con->compute_normal();
//    }
}


void Mesh::build_connection(int i, int j, int k, int h, \
                            int opposite, Element* elem,  \
                            std::map<std::array<unsigned int,4>, Connection*> &unique_id_map,
                            std::vector<Connection*> &temp_connection) {

    std::array<unsigned int, 4> key = {0,0,0,0};
    std::map<std::array<unsigned int, 4>, Connection*>::iterator it;
    Connection* con;

    //build key
    key[0] = elem->vertice_ids[i];
    key[1] = elem->vertice_ids[j];
#ifdef DEBUG_MODE
    for (auto x : key) std::cout << x << ' ';
    std::cout << std::endl;
#endif
    if (k>0) key[2] = elem->vertice_ids[k]; //if 3D mesh
    if (h>0) key[3] = elem->vertice_ids[h]; 
    std::sort(key.begin(),key.end());
#ifdef DEBUG_MODE
    for (auto x : key) std::cout << x << ' ';
    std::cout << std::endl << std::endl;
#endif

    //find connection in map
    it = unique_id_map.find(key);
    if (it != unique_id_map.end()) { //key already in
        //terminate connection connection
        it->second->element_id_up = elem;
        it->second->vertice_up = elem->vertices[opposite];
        it->second->check_orientation();
        unique_id_map.erase(it);
    }
    else {
        con = new Connection();
        con->vertices.push_back(elem->vertices[i]);
        con->vertices.push_back(elem->vertices[j]);
        if (dim == 3) con->vertices.push_back(elem->vertices[k]);
        if (h>0) con->vertices.push_back(elem->vertices[h]);
        con->element_id_dn = elem;
        con->vertice_dn = elem->vertices[opposite];
        unique_id_map[key] = con;
        temp_connection.push_back(con);
    }
}


// Load from array
void Mesh::load_vertices_array(double vertices[], size_t n_vertices) {
    //load vertices from a array
    //designed for Python external call
    //vertices array must be 1D of size n_vertices*3
    for (size_t id=0; id!=n_vertices; id++) {
        add_vertice(vertices[3*id], vertices[3*id+1],
                          vertices[3*id+2], id+1);
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
            if (elements[count] <= 0) {continue; }
            ids.push_back(elements[count]);
            count++;
        }
        add_element(ids); //create the element
    }
}


void Mesh::save_face_non_orthogonality_angle(std::string f_out) {
    std::ofstream out(f_out);
    unsigned int inverted_element = 0;
    for (Connection* con : connections_internal) {
        out << con->element_id_dn->natural_id << " ";
        out << con->element_id_up->natural_id << " ";
        //out << std::acos(std::abs(con->orthogonality)) * 57.29583 << std::endl; //57.2583 = 180 / pi
        out << 1-con->compute_orthogonality() << std::endl;
    	if (con->orthogonality < 0) {
    	    inverted_element++;
    	}
    }
    if (inverted_element != 0) {
        std::cout << "WARNING: " <<  inverted_element;
        std::cout << " inverted element found!" << std::endl;
    }
    out.close();
}

void Mesh::save_face_detailed_informations(std::string f_out) {
    std::ofstream out(f_out);
    out << "id1 id2 v1 v2 v3 v4 area nx ny nz cx cy cz error" << std::endl;
    for (Connection* con : connections_internal) {
        out << con->element_id_dn->natural_id << " ";
        out << con->element_id_up->natural_id << " ";
        for (Vertice* v: con->vertices) {out << v->natural_id << " ";}
        if (con->vertices.size() == 3) {out << "-1 ";}
        out << con->area << " ";
        for (size_t i=0; i<dim; i++) out << (con->normal)[i] << " ";
        for (size_t i=0; i<dim; i++) out << (con->cell_center_vector)[i] << " ";
        con->compute_orthogonality();
        
        out << 1-con->orthogonality << " " << std::endl;
    }
    out.close();
}

void Mesh::display_stats() {
    double mean = 0;
    double max = -999;
    double error;
    for (Connection* con : connections_internal) {
        error = 1-con->orthogonality;
        if (error > max) {max = error;};
        mean += error;
    }
    mean /= n_connections_internal;
    //maximum
    std::cout << "Maximum error: " << max << std::endl;
    //mean
    std::cout << "Mean error: " << mean << std::endl;
}
