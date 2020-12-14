#ifndef BINDING_H_INCLUDED
#define BINDING_H_INCLUDED

#include "../cpp/OrthOpt.h"
#include "../cpp/Mesh.h"
#include <iostream>

extern "C" {
    //LOAD MESH
    Mesh* Mesh_init() {return new Mesh(); }
    void Mesh_load_vertices_array(Mesh* mesh, double* vertices[], size_t s) {
        std::cout << *vertices[0] << std::endl;
        mesh->load_vertices_array(*vertices, s);
    }
    void Mesh_load_elements_array(Mesh* mesh, unsigned int* elements[],
                                  unsigned int* types[], size_t s) {
        mesh->load_elements_array(*elements, *types, s);
    }
    void Mesh_add_vertice(Mesh* mesh, double x, double y, 
                          double z, unsigned int id) {
        mesh->add_vertice(x,y,z,id);
    }
    void Mesh_add_element(Mesh* mesh, unsigned int* ids[]) {
        //mesh->add_element(*ids);
    }
    
    //CREATE ORTHOPT OBJECT
    OrthOpt* OrthOpt_init() {return new OrthOpt(); };
    void OrthOpt_set_penalizing_power(OrthOpt* obj, int p) {
        obj->set_penalizing_power(p);
    }
    unsigned int optimize(OrthOpt* obj, double penalization_power,
                          unsigned int max_it) {
        unsigned int code = 0;
        
        return code;
    }
    
    int test() {return 1;}
    
}
 



#endif // BINDING_H_INCLUDED
