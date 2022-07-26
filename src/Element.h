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

        Element() {};
        virtual ~Element() {};

        virtual Point center();
        
        virtual Eigen::Matrix3d center_derivative(Vertice* p);
};

class Tri_Element : public Element
{
    public:
        Tri_Element() {
            type = -3;
        }
        Point Center() {
            Point center(0., 0., 0.);
            for (auto p = vertices.begin(); p != vertices.end(); p++) {
                center += *((*p)->coor);
            }
            center /= 3;
            return center;
        }
        
        Eigen::Matrix3d center_derivative(Vertice* p) {
            Eigen::Matrix3d derivative;
            derivative << 0.333333, 0, 0, 0, 0.333333, 0, 0, 0, 0;
            return derivative;
        }
};

class Tet_Element : public Element
{
    public:
        Tet_Element() {
            type = 4;
        }
        
        Point Center() {
            Point center(0., 0., 0.);
            for (auto p = vertices.begin(); p != vertices.end(); p++) {
                center += *((*p)->coor);
            }
            center /= 4;
            return center;
        }
        
        Eigen::Matrix3d center_derivative(Vertice* p) {
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
        
        Point Center() {
            Point center(0., 0., 0.);
            //WARNING, polygon vertices must be ordered and define the polygon segment
            Point* vi = nullptr; 
            Point *vip = nullptr;
            double temp, area = 0.;
            for (size_t i=0; i!=vertices.size(); i++) {
                vi = vertices[i]->coor; 
                if (i+1 == vertices.size()) vip = vertices[0]->coor;
                else vip = vertices[i+1]->coor;
                temp = (vi->x*vip->y - vip->x*vi->y);
                center.x += (vi->x + vip->x)*temp;
                center.y += (vi->y + vip->y)*temp;
                area += temp;
            }
            center /= 6*area;
            return center;
        }
        
        Eigen::Matrix3d center_derivative(Vertice* p) {
            Eigen::Matrix3d derivative;
            //calculate center and area
            Point center(0., 0., 0.);
            Point* vi = nullptr; 
            Point *vip = nullptr;
            double temp, area = 0.;
            for (size_t i=0; i!=vertices.size(); i++) {
                vi = vertices[i]->coor; 
                if (i+1 == vertices.size()) vip = vertices[0]->coor;
                else vip = vertices[i+1]->coor;
                temp = (vi->x*vip->y - vip->x*vi->y);
                center.x += (vi->x + vip->x)*temp;
                center.y += (vi->y + vip->y)*temp;
                area += temp;
            }
            center /= 6*area;
            //derivative
            double dxx, dxy, dyx, dyy;
            double xm, ym, x, y, xp, yp;
            dxx = (xm * (y - ym) + 2*x*(yp-ym) + xp*(yp-y)) / 6 + center.x * 0.5 * (yp-ym);; 
            dxy = (xm*xm + x*(xm-xp) + xp*xp) / 6 + center.x * 0.5 * (xp-xm); 
            dyx = (-ym*ym + y*(yp-ym) + yp*yp) / 6 + center.y * 0.5 * (yp-ym); 
            dxx = (ym * (x - xp) + 2*y*(xm-xp) + yp*(xm-x)) / 6 + center.y * 0.5 * (xp-xm); 
            derivative << dxx, dxy, 0, dyx, dyy, 0, 0, 0, 0;
            derivative /= area;
            return derivative;
        }
};



#endif // ELEMENT_H
