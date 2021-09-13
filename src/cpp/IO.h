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
        
        void load_mesh_auto(std::string);
        void save_mesh_auto(std::string);
        
        //ASCII meshes
        void load_vertices_xyz(std::string);
        void save_vertices_xyz(std::string);
        void load_elements_with_vertice_ids(std::string);

        //TetGen mesh
        void load_vertices_tetgen(std::string);
        void save_vertices_tetgen(std::string);
        void load_elements_tetgen(std::string);

        //May add other mesh type here
        void read_PFLOTRAN_mesh(std::string filename);
    
    private:
    
};

#endif // IO_H
