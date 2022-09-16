#ifndef VERTICE_H
#define VERTICE_H

#include <Eigen/Core>
#include <iostream>

class Vertice
{
    public:
        Eigen::VectorXd *coor;
        unsigned int natural_id=-1;
        bool fixed = false;
        
        Vertice(double x, double y) {
            coor = new Eigen::VectorXd(2);
            (*coor)[0] = x;
            (*coor)[1] = y;
        };
        Vertice(double x, double y, unsigned int id) {
            coor = new Eigen::VectorXd(2);
            (*coor)[0] = x;
            (*coor)[1] = y;
            natural_id = id;
        };
        Vertice(double x, double y, double z) {
            coor = new Eigen::VectorXd(3);
            (*coor)[0] = x;
            (*coor)[1] = y;
            (*coor)[2] = z;
        };
        Vertice(double x, double y, double z, unsigned int id) {
            coor = new Eigen::VectorXd(3);
            (*coor)[0] = x;
            (*coor)[1] = y;
            (*coor)[2] = z;
            natural_id = id;
        };
        virtual ~Vertice () {};
        void fix_vertice() {fixed = true;}
        
        friend std::ostream& operator<<(std::ostream& os, const Vertice& v) {
            //to print connection information (debug purpose)
            os << "vertices id: " << v.natural_id << " ; ";
            os << "fixed: " << v.fixed << " ; " ;
            os << "coor: " << v.coor;
            return os;
        }
};

#endif // VERTICE_H
