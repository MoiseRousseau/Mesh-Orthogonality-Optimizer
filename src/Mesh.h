#ifndef MESH_H
#define MESH_H

#include <string>
#include <map>
#include <array>
#include "Connection.h"
#include "Element.h"
#include "Vertice.h"


class Mesh
{
    public:

        unsigned int dim = 3;
        unsigned int n_vertices = 0;
        std::vector<Vertice*> vertices;
        unsigned int n_elements = 0;
        std::vector<Element*> elements;
        
        unsigned int n_connections_internal = 0;
        std::vector<Connection*> connections_internal;
        unsigned int n_connections_bc = 0;
        std::vector<Connection*> boundary_connections;

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
            for (auto con = connections_internal.begin(); 
                con != connections_internal.end(); con++) {
                delete *con;
            }
            connections_internal.clear();
            for (auto con = boundary_connections.begin(); 
                con != boundary_connections.end(); con++) {
                delete *con;
            }
            boundary_connections.clear();
        };
        
        //adding one vertices and one element at a time
        void add_vertice(double x, double y, unsigned int id) {
            //note here: we explicitely ask for vertice id since they are used
            // to defined the elements.
            n_vertices += 1;
            Vertice* A = new Vertice(x,y,id);
            vertices.push_back(A);
        }
        void add_vertice(double x, double y, double z, unsigned int id) {
            //note here: we explicitely ask for vertice id since they are used
            // to defined the elements.
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
            if (dim == 2) elem->type = -ids.size();
            else elem->type = ids.size();
            elem->natural_id = n_elements;
            elements.push_back(elem);
        }
        
        
        //Add binding with python and numpy
        void load_vertices_array(double vertices[], size_t n_vertices);
        void load_elements_array(unsigned int elements[], unsigned int type[],
                                 size_t s);
        double* save_elements_array();
        
        
        /**
         * \brief Decompose the mesh to find internal and boundary
              connections
         * \details The part of the input mesh should be specified as
         *    a contiguous range of facet indices.
         * \param[in] facets_begin first facet in the range
         * \param[in] facets_end one past last facet in the range
         */
        void decompose();
        
        
        //output
        void save_face_non_orthogonality_angle(std::string f_out);
        void save_face_detailed_informations(std::string f_out);
        void display_stats();


    private:
        void build_connection(int i, int j, \
                              int k, int h, \
                              int opposite, Element* elem,  \
                              std::map<std::array<unsigned int, 4>, Connection*> &unique_id_map,
                              std::vector<Connection*> &temp);
};

#endif // MESH_H
