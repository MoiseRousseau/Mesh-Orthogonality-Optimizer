#ifndef VERTICE_H
#define VERTICE_H

#include <Eigen/Core>

class Vertice
{
    public:
        Eigen::VectorXd* coor;
        unsigned int natural_id=-1;
        bool fixed = false;
        
        Vertice(double x, double y) {
            coor = new Eigen::VectorXd;
            coor->resize(2);
            *coor << x,y;
        };
        Vertice(double x, double y, unsigned int id) {
            Vertice(x,y);
            natural_id = id;
        };
        Vertice(double x, double y, double z) {
            coor = new Eigen::VectorXd;
            coor->resize(3);
            *coor << x,y,z;
        };
        Vertice(double x, double y, double z, unsigned int id) {
            Vertice(x,y,z);
            natural_id = id;
        };
        virtual ~Vertice () {};
        void fix_vertice() {fixed = true;}
};

#endif // VERTICE_H
