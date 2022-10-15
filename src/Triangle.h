#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"
#include <cfloat>

struct RayTriangleIntersection{
    bool intersectionExists;
    float t;
    float w0,w1,w2; // Coordonnées barycentriques de l'intersection avec le triangle
    unsigned int tIndex;
    Vec3 intersection;
    Vec3 normal;
    RayTriangleIntersection() : intersectionExists(false) , t(FLT_MAX) {}
};

class Triangle {
private:
    Vec3 m_c[3] , m_normal;
    float area;
public:
    Triangle() {}
    Triangle( Vec3 const & c0 , Vec3 const & c1 , Vec3 const & c2 ) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
    }
    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross( m_c[1] - m_c[0] , m_c[2] - m_c[0] );
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }
    void setC0( Vec3 const & c0 ) { m_c[0] = c0; } // remember to update the area and normal afterwards!
    void setC1( Vec3 const & c1 ) { m_c[1] = c1; } // remember to update the area and normal afterwards!
    void setC2( Vec3 const & c2 ) { m_c[2] = c2; } // remember to update the area and normal afterwards!
    Vec3 const & normal() const { return m_normal; }
    Vec3 projectOnSupportPlane( Vec3 const & p ) const {
        Vec3 result;
        //TODO completer
        return result;
    }
    float squareDistanceToSupportPlane( Vec3 const & p ) const {
        float result;
        //TODO completer
        return result;
    }
    float distanceToSupportPlane( Vec3 const & p ) const { return sqrt( squareDistanceToSupportPlane(p) ); }
    bool isParallelTo( Line const & L ) const {
        //TODO completer
        return Vec3::dot(m_normal, L.direction()) == 0;
    }
    Vec3 getIntersectionPointWithSupportPlane( Line const & L ) const {
        // you should check first that the line is not parallel to the plane!
        Vec3 result;
        //TODO completer
        return result;
    }

    RayTriangleIntersection getIntersection( Ray const & ray ) const {

        Vec3 edge1;
        Vec3 edge2;
        Vec3 n;
        float t;
        Vec3 rayPoint;
        RayTriangleIntersection result;

        if(!isParallelTo(ray)) {

            edge1 = this->m_c[2] - this->m_c[0];
            edge2 = this->m_c[1] - this->m_c[0];
            n = Vec3::cross(edge2,edge1);
            n.normalize();

            t = (Vec3::dot(this->m_c[0],n) - Vec3::dot(ray.origin(),n)) / Vec3::dot(ray.direction(),n);

            if(t > 0.00001) {

                rayPoint = ray.origin() + t * ray.direction();

                result.w0 = Vec3::dot(this->m_normal,Vec3::cross(edge2,rayPoint - this->m_c[0]));
                result.w1 = Vec3::dot(this->m_normal,Vec3::cross(this->m_c[2] - this->m_c[1],rayPoint - this->m_c[1]));
                result.w2 = Vec3::dot(this->m_normal,Vec3::cross(this->m_c[0] - this->m_c[2],rayPoint - this->m_c[2]));

                if(result.w0 >= 0 && result.w1 >= 0 && result.w2 >= 0)  {

                    result.intersectionExists = true;
                    result.intersection = rayPoint;
                    result.normal = n;
                    result.t = t;
                }
            }
        }

        return result;
    }
};
#endif
