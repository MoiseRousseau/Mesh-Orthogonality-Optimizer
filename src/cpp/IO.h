#ifndef IO_H
#define IO_H

#include "Mesh.h"
#include <string>


class IO
{
    public:
    
        Mesh* mesh; //the mesh for which IO is directed
        
        IO(Mesh* m) {
            mesh = m;
        }
        
        void load_mesh_auto(std::string, std::string, std::string);
        void save_mesh_auto(std::string);
        
        //ASCII meshes
        void load_vertices_xyz(std::string);
        void save_vertices_xyz(std::string);
        void load_elements_with_vertice_ids(std::string);

        //TetGen mesh
        void load_vertices_tetgen(std::string);
        void save_vertices_tetgen(std::string);
        void load_elements_tetgen(std::string);
        
        //Medit mesh
        void load_mesh_medit(std::string);
        void save_mesh_medit(std::string);
        
        //PFLOTRAN mesh
        void load_mesh_PFLOTRAN(std::string);
        void save_mesh_PFLOTRAN(std::string);
    
    private:
        //Medit utils
        bool s_eqi(std::string, std::string);
        bool s_begin(std::string, std::string);
        int s_len_trim(std::string);
};

#endif // IO_H
