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
    std::map<std::array<unsigned int, 4>, Connection*> unique_id_map;
    
    //link element to vertice instances
    for (Element* elem : elements) {
        for (unsigned int id : elem->vertice_ids) {
            elem->vertices.push_back(vertices[id-1]);
        }
    }
    for (Element* elem : elements) {
        if (elem->type == 4) { //tetrahedron
            //first face 0 1 2, opposite 3
            build_connection(0,1,2,-1, 3, elem, unique_id_map);

            //first face 1 2 3, opposite 0
            build_connection(1,2,3,-1,0, elem, unique_id_map);

            //first face 0 1 3, opposite 2
            build_connection(0,1,3,-1, 2, elem, unique_id_map);

            //first face 0 2 3, opposite 1
            build_connection(0,2,3,-1,1, elem, unique_id_map);
        }
        else if (elem->type == 5) { //pyramid
            build_connection(0,1,2,3,4, elem, unique_id_map);
            //other are treated by tet above
        }
        else if (elem->type == 6) { //prisms
            //not optimized, so do not figure
        }
        else if (elem->type == 8) { //hex
            //not optimized, so do not figure
        }
        else {
            exit(1);
        }
    }
    for (auto con = connections_internal.begin(); 
         con != connections_internal.end();) {
        if ((*con)->vertice_up == nullptr) {
            for (Vertice* v : (*con)->vertices) {
                v->fix_vertice();
            }
            connections_internal.erase(con);
            boundary_connections.push_back(*con);
        }
        else {++con;}
    }
    n_connections_internal = connections_internal.size();
    n_connections_bc = boundary_connections.size();
    
    for (Connection* con : connections_internal) {
        con->compute_normal();
    }
}


void Mesh::build_connection(int i, int j, int k, int h, \
                            int opposite, Element* elem,  \
                            std::map<std::array<unsigned int, 4>, Connection*> &unique_id_map) {

    std::array<unsigned int, 4> key;
    std::map<std::array<unsigned int, 4>, Connection*>::iterator it;
    Connection* con;

    //build key
    key[0] = elem->vertice_ids[i];
    key[1] = elem->vertice_ids[j];
    key[2] = elem->vertice_ids[k];
    if (h>0) {key[3] = elem->vertice_ids[h];}
    else {key[3] = 0;}
    std::sort(key.begin(),key.end());

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
        con->vertices.push_back(elem->vertices[k]);
        con->type = 3;
        if (h>0) {con->vertices.push_back(elem->vertices[h]); con->type++;}
        con->element_id_dn = elem;
        con->vertice_dn = elem->vertices[opposite];
        unique_id_map[key] = con;
        connections_internal.push_back(con);
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
    for (Connection* con : connections_internal) {
        out << con->element_id_dn->natural_id << " ";
        out << con->element_id_up->natural_id << " ";
        //out << std::acos(std::abs(1-con->error)) * 57.29583 << std::endl; //57.2583 = 180 / pi
        out << con->error << std::endl;
    }
    out.close();
}

void Mesh::save_face_informations(std::string f_out) {
    std::ofstream out(f_out);
    out << "id1 id2 v1 v2 v3 v4 area nx ny nz cx cy cz error\n";
    for (Connection* con : connections_internal) {
        out << con->element_id_dn->natural_id << " ";
        out << con->element_id_up->natural_id << " ";
        for (Vertice* v: con->vertices) {out << v->natural_id << " ";}
        if (con->vertices.size() == 3) {out << "-1 ";}
        out << con->area << " ";
        out << con->normal.x << " " << con->normal.y << " " << con->normal.z << " ";
        out << con->cell_center_vector.x << " " << con->cell_center_vector.y << " " << con->cell_center_vector.z << " ";
        out << con->error << " " << std::endl;
    }
    out.close();
}

void Mesh::display_stats() {
    double mean = 0;
    double max = -999;
    for (Connection* con : connections_internal) {
        if (con->error > max) {max = con->error;};
        mean += con->error;
    }
    mean /= n_connections_internal;
    //maximum
    std::cout << "Maximum error: " << max << std::endl;
    //mean
    std::cout << "Mean error: " << mean << std::endl;
}
