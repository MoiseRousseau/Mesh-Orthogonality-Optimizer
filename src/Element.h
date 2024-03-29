#ifndef ELEMENT_H
#define ELEMENT_H

#include "Vertice.h"
#include <Eigen/Core>


class Element
{
    public:

        unsigned int natural_id=0;
        std::vector<unsigned int> vertice_ids; //TODO: not needed, already in Vertices*
        int type; //4 = tet, 5 = pyr, 6 = prisms, 8 = hex, -3 triangle, -4 = quad
        std::vector<Vertice*> vertices;
        int zone = 0; //material

        Element() {};
        virtual ~Element() {};
        
        virtual Eigen::VectorXd center() {return Eigen::VectorXd();};
        //virtual Eigen::Vector3d center();
        
        virtual Eigen::MatrixXd center_derivative(Vertice* p) {return Eigen::MatrixXd();};
};

class Tri_Element : public Element
{
    public:
        Tri_Element() {
            type = -3;
        }
        Eigen::VectorXd center() {
            Eigen::Vector2d center(0., 0.);
            for (auto p = vertices.begin(); p != vertices.end(); p++) {
                center += *((*p)->coor);
            }
            center /= 3;
            return center;
        }
        
        Eigen::MatrixXd center_derivative(Vertice* p) {
            Eigen::Matrix2d derivative;
            derivative << 0.333333, 0, 0, 0.333333;
            return derivative;
        }
};

class Tet_Element : public Element
{
    public:
        Tet_Element() {
            type = 4;
        }
        
        Eigen::VectorXd center() {
            Eigen::Vector3d center(0., 0., 0.);
            for (auto p = vertices.begin(); p != vertices.end(); p++) {
                center += *((*p)->coor);
            }
            center /= 4;
            return center;
        }
        
        Eigen::MatrixXd center_derivative(Vertice* p) {
            Eigen::Matrix3d derivative;
            derivative << 0.25, 0, 0, 0, 0.25, 0, 0, 0, 0.25;
            return derivative;
        }
};

class Polygon_Element : public Element
{
    public:
    
        Polygon_Element() {
            type = -4;
        }
        
        Eigen::VectorXd center() {
            Eigen::Vector2d center(0., 0.);
            //WARNING, polygon vertices must be ordered
            Eigen::VectorXd *vi = nullptr; 
            Eigen::VectorXd *vip = nullptr;
            double temp, area = 0.;
            for (size_t i=0; i!=vertices.size(); i++) {
                vi = vertices[i]->coor; 
                if (i+1 == vertices.size()) vip = vertices[0]->coor;
                else vip = vertices[i+1]->coor;
                temp = (*vi)[0] * (*vip)[1] - (*vip)[0] * (*vi)[1] ;
                center[0] += ( (*vi)[0] + (*vip)[0]) * temp;
                center[1] += ((*vi)[1] + (*vip)[1])*temp;
                area += temp/2;
            }
            center /= 6*area;
            return center;
        }
       
        Eigen::MatrixXd center_derivative(Vertice* p) {
#if 0 
            Eigen::Matrix2d derivative;
            //calculate center and area
            Eigen::Vector2d center(0., 0.);
            Eigen::VectorXd *vi = nullptr;
            Eigen::VectorXd *vip = nullptr;
            double temp, area = 0.;
            int index=0;
            for (size_t i=0; i!=vertices.size(); i++) {
                if (vertices[i] == p) index = i;
                vi = vertices[i]->coor; 
                if (i+1 == vertices.size()) vip = vertices[0]->coor;
                else vip = vertices[i+1]->coor;
                temp = ((*vi)[0] * (*vip)[1] - (*vip)[0] * (*vi)[1]);
                //center[0] += (vi[0] + vip[0])*temp;
                //center[1] += (vi[1] + vip[1])*temp;
                center += (*vi + *vip) * temp;
                area += temp/2;
            }
            center /= 6*area;
            
            //derivative
            double dxx, dxy, dyx, dyy;
            double xm, ym, x, y, xp, yp;
            //vi
            x = (*vertices[index]->coor)[0];
            y = (*vertices[index]->coor)[1];
            //vim
            if (index == 0) {
                xm = (*vertices.back()->coor)[0];
                ym = (*vertices.back()->coor)[1];
            }
            else {
                xm = (*vertices[index-1]->coor)[0];
                ym = (*vertices[index-1]->coor)[1];
            }
            //vip
            if (index+1 == vertices.size()) {
                xp = (*vertices[0]->coor)[0];
                yp = (*vertices[0]->coor)[1];
            }
            else {
                xp = (*vertices[index+1]->coor)[0];
                yp = (*vertices[index+1]->coor)[1];
            }
            
            dxx = (xm * (y - ym) + 2*x*(yp-ym) + xp*(yp-y)) / 6 + center[0] * 0.5 * (yp-ym); 
            dxy = (xm*xm + x*(xm-xp) + xp*xp) / 6 + center[0] * 0.5 * (xp-xm); 
            dyx = (-ym*ym + y*(yp-ym) + yp*yp) / 6 + center[1] * 0.5 * (yp-ym); 
            dxx = (ym * (x - xp) + 2*y*(xm-xp) + yp*(xm-x)) / 6 + center[1] * 0.5 * (xp-xm); 
            derivative << dxx, dxy, dyx, dyy;
            derivative /= area;
            std::cout << "deriv" << std::endl << derivative << std::endl << "end" << std::endl;
            return derivative;
#endif
            //TODO: this is not thread safe!!!
            double pertub = 1e-8;
            double old_coor;

            Eigen::Vector2d c = center();
            
            old_coor = (*p->coor)[0];
            (*p->coor)[0] += pertub;
            Eigen::Vector2d c_dx = center();
            (*p->coor)[0] = old_coor;
            
            old_coor = (*p->coor)[1];
            (*p->coor)[1] += pertub;
            Eigen::Vector2d c_dy = center();
            (*p->coor)[1] = old_coor;
            
            Eigen::Matrix2d derivative;
            derivative << (c_dx[0]-c[0]) / pertub, (c_dy[0]-c[0]) / pertub, (c_dx[1]-c[1]) / pertub, (c_dy[1]-c[1]) / pertub;
            return derivative;
        }
};



#endif // ELEMENT_H
