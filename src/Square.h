#ifndef SQUARE_H
#define SQUARE_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySquareIntersection{
    bool intersectionExists;
    float t;
    float u,v;
    Vec3 intersection;
    Vec3 normal;
    RaySquareIntersection() : intersectionExists(false), t(FLT_MAX) {}
};


class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    Square() : Mesh() {}
    Square(Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
           float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void setQuad( Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
                  float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector , upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector*width;
        m_up_vector = m_up_vector*height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;                                      vertices[0].u = uMin; vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;                     vertices[1].u = uMax; vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;       vertices[2].u = uMax; vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;                        vertices[3].u = uMin; vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;
    }

    RaySquareIntersection intersect(const Ray &ray) const {
        
        //TODO calculer l'intersection rayon quad

        float D;
        float t;
        Vec3 a;
        Vec3 n;
        Vec3 p;
        Vec3 point1;
        Vec3 point2;
        Vec3 point3;
        Vec3 point4;
        RaySquareIntersection raySquareIntersection;

        // D'après la formule du cours de l'intersection rayon-plan page 44, on commence par déterminer les coefficients calculables pour exprimer en fonction de t
        a = this->vertices[0].position;
        n = Vec3::cross(this->vertices[1].position - a,this->vertices[3].position - a);
        n.normalize();
        D = Vec3::dot(a,n);

        t = (D - Vec3::dot(ray.origin(),n)) / Vec3::dot(ray.direction(),n);

        p = ray.origin() + t * ray.direction();

        // Puis on récupère l'intersection seulement si t est positif
        if(t > 0) {

            raySquareIntersection.intersectionExists = true;

            point1 = Vec3::cross(this->vertices[1].position - this->vertices[0].position,p - this->vertices[0].position);
            point2 = Vec3::cross(this->vertices[2].position - this->vertices[1].position,p - this->vertices[1].position);
            point3 = Vec3::cross(this->vertices[3].position - this->vertices[2].position,p - this->vertices[2].position);
            point4 = Vec3::cross(this->vertices[0].position - this->vertices[3].position,p - this->vertices[3].position);

            if( ((Vec3::dot(point1,n) > 0) == (Vec3::dot(point2,n) > 0)) && ((Vec3::dot(point1,n) > 0) == (Vec3::dot(point3,n) > 0)) && ((Vec3::dot(point1,n) > 0) == (Vec3::dot(point4,n) > 0)) ) {

                raySquareIntersection.t = t;
                raySquareIntersection.intersection = p; 
                raySquareIntersection.normal = n;
            }
        }

        return raySquareIntersection;
    }
};
#endif // SQUARE_H
