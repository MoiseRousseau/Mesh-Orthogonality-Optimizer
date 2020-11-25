#ifndef MESH_H
#define MESH_H

#include <string>
#include "Connection.h"
#include "Point.h"
#include "Element.h"
#include "Vertice.h"


class Mesh
{
    public:

        std::vector<Vertice*> vertices;
        std::vector<Element*> elements;
        unsigned int n_vertices = 0;
        unsigned int n_elements = 0;

        Mesh() {};

        virtual ~Mesh() {
            for (auto v = vertices.begin(); v != vertices.end(); v++) {
                delete *v;
            }
            vertices.clear();
            for (auto v = elements.begin(); v != elements.end(); v++) {
                delete *v;
            }
            elements.clear();
        };

        // // // INPUT // // //
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

        //Add binding with python and numpy
        void read_numpy_vertices() {};
        void read_numpy_elements() {};

        //adding one vertices and one element at a time
        void add_vertice(double x, double y, double z, unsigned int id) {
            n_vertices += 1;
            Vertice* A = new Vertice(x,y,z,id);
            vertices.push_back(A);
        }
        void add_element(std::vector<unsigned int> ids) {
            n_elements += 1;
            Element* elem = new Element();
            for (auto id = ids.begin(); id != ids.end(); id++) {
                elem->vertice_ids.push_back(*id);
            }
            elem->type = ids.size();
            elem->natural_id = n_elements;
            elements.push_back(elem);
        }

        void decompose_mesh() {
            //link element to vertice instances
            for (Element* elem : elements) {
                for (unsigned int id : elem->vertice_ids) {
                    elem->vertices.push_back(vertices[id-1]);
                }
            }
            //TODO delete unused vertices ?
        }

    private:
};

#endif // MESH_H
