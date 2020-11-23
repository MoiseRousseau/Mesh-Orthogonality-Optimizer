#ifndef ELEMENT_H
#define ELEMENT_H

#include "Vertice.h"


class Element
{
    public:

        unsigned int natural_id;
        std::vector<unsigned int> vertice_ids; //from the input
        unsigned int type; //4 = tet, 5 = pyr, 6 = prisms, 8 = hex
        std::vector<Vertice*> vertices;

        Element() {
            natural_id = 0;
            type = 0;
        };
        virtual ~Element() {};

        Point center() {
            Point center(0., 0., 0.);
            for (auto p = vertices.begin(); p != vertices.end(); p++) {
                center += *((*p)->coor);
            }
            center /= type; //type is the number of vertice
            return center;
        }
};



#endif // ELEMENT_H
