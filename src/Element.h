#ifndef ELEMENT_H
#define ELEMENT_H

#include "Vertice.h"
#include <Eigen/Core>


class Element
{
    public:

        unsigned int natural_id;
        std::vector<unsigned int> vertice_ids; //from the input
        int type; //4 = tet, 5 = pyr, 6 = prisms, 8 = hex, -3 triangle, -4 = quad
        std::vector<Vertice*> vertices;

        Element() {
            natural_id = 0;
            type = 0;
        };
        virtual ~Element() {};

        Eigen::Vector3d center() {
            Eigen::Vector3d center(0., 0., 0.);
            if (type == -3 or type == 4) { 
                //tet or triangle, center = mean of coordinate
                for (auto p = vertices.begin(); p != vertices.end(); p++) {
                    center += *((*p)->coor);
                }
                center /= vertice_ids.size();
            }
            else if (type <= -4) { //other polygon including quad
                //WARNING, polygon vertices must be ordered and define the polygon segment
                Eigen::Vector3d* vi = nullptr; 
                Eigen::Vector3d *vip = nullptr;
                double temp, area = 0.;
                for (size_t i=0; i!=vertices.size(); i++) {
                    vi = vertices[i]->coor; 
                    if (i+1 == vertices.size()) vip = vertices[0]->coor;
                    else vip = vertices[i+1]->coor;
                    temp = (*vi)[0] * (*vip)[1] - (*vip)[0] * (*vi)[1] ;
                    center[0] += ( (*vi)[0] + (*vip)[0]) * temp;
                    center[1] += ((*vi)[1] + (*vip)[1])*temp;
                    area += temp;
                }
                center /= 6*area;
            }
            return center;
        }
};



#endif // ELEMENT_H
