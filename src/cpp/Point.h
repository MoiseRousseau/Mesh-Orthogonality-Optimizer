#ifndef POINT_H
#define POINT_H
#include <cmath>
#include <iostream>

class Point
{
    public:

        double x, y, z;

        Point() {x = 0.; y = 0.; z = 0;};
        Point(double a, double b, double c) {
            x = a; y = b; z = c;
        }
        virtual ~Point() {};

        Point operator+ (Point X) {return Point(x + X.x, y + X.y, z + X.z);}
        Point operator+ (double a) {return Point(a+x, a+y, a+z);}
        Point operator/ (double n) {return Point(x/n, y/n, z/n);}
        Point operator- (Point X) {return Point(x-X.x, y-X.y, z-X.z);}
        Point operator- (void) {return Point(-x, -y, -z);}
        Point operator* (double a) {return Point(a*x, a*y, a*z);}

        //Point operator= (double n) {return Point(n, n, n);}
        void operator= (double a) {x = a; y = a; z = a;}
        //Point operator= (Point A) {return Point(A.x, A.y, A.z);}
        void operator+= (Point A) {x += A.x; y += A.y; z += A.z;}
        void operator/= (double n) {x /= n; y /= n; z /= n;}

        double dot(Point A) {return x*A.x + y*A.y + z*A.z;}
        double norm() {return std::pow(x*x+y*y+z*z,0.5);}
        Point cross(Point A) {
            return Point(y*A.z - A.y*z,
                         z*A.x - A.z*x,
                         x*A.y - A.x*y);
        }
        Point cross(Point* A) {
            return Point(y*A->z - A->y*z,
                         z*A->x - A->z*x,
                         x*A->y - A->x*y);
        }

    protected:

    private:
};

#endif // POINT_H
