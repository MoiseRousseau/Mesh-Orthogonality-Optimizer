#ifndef VERTICE_H
#define VERTICE_H

#include <Eigen/Core>

class Vertice
{
    public:
        Eigen::Vector3d* coor;
        unsigned int natural_id;
        bool fixed = false;
        
        Vertice(double x, double y) {
            coor = new Eigen::Vector3d(x,y,0.);
            natural_id = -1;
        };
        Vertice(double x, double y, unsigned int id) {
            coor = new Eigen::Vector3d(x,y,0.);
            natural_id = id;
        };
        Vertice(double x, double y, double z) {
            coor = new Eigen::Vector3d(x,y,z);
            natural_id = -1;
        };
        Vertice(double x, double y, double z, unsigned int id) {
            coor = new Eigen::Vector3d(x,y,z);
            natural_id = id;
        };
        virtual ~Vertice () {};
        void fix_vertice() {fixed = true;}
};

#endif // VERTICE_H
