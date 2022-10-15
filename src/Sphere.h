#ifndef Sphere_H
#define Sphere_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySphereIntersection{
    bool intersectionExists;
    float t;
    float theta,phi; // Coordonnées sphériques du point d'intersection, pour le plaquage de textures plus loin dans le projet
    unsigned int sIndex;
    Vec3 intersection;
    Vec3 secondintersection;
    Vec3 normal; // Calcul du vecteur du point d'intersection au centre de la sphère et normaliser
    RaySphereIntersection() : intersectionExists(false) , t(FLT_MAX) {}
};

static
Vec3 SphericalCoordinatesToEuclidean( Vec3 ThetaPhiR ) {
    return ThetaPhiR[2] * Vec3( cos(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[1]) );
}
static
Vec3 SphericalCoordinatesToEuclidean( float theta , float phi ) {
    return Vec3( cos(theta) * cos(phi) , sin(theta) * cos(phi) , sin(phi) );
}

static
Vec3 EuclideanCoordinatesToSpherical( Vec3 xyz ) {
    float R = xyz.length();
    float phi = asin( xyz[2] / R );
    float theta = atan2( xyz[1] , xyz[0] );
    return Vec3( theta , phi , R );
}



class Sphere : public Mesh {
public:
    Vec3 m_center;
    float m_radius;

    Sphere() : Mesh() {}
    Sphere(Vec3 c , float r) : Mesh() , m_center(c) , m_radius(r) {}

    void build_arrays(){
        unsigned int nTheta = 20 , nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi );
        normalsArray.resize(3 * nTheta * nPhi );
        uvs_array.resize(2 * nTheta * nPhi );
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta ; ++thetaIt ) {
            float u = (float)(thetaIt) / (float)(nTheta-1);
            float theta = u * 2 * M_PI;
            for( unsigned int phiIt = 0 ; phiIt < nPhi ; ++phiIt ) {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float)(phiIt) / (float)(nPhi-1);
                float phi = - M_PI/2.0 + v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean( theta , phi );
                positions_array[ 3 * vertexIndex + 0 ] = m_center[0] + m_radius * xyz[0];
                positions_array[ 3 * vertexIndex + 1 ] = m_center[1] + m_radius * xyz[1];
                positions_array[ 3 * vertexIndex + 2 ] = m_center[2] + m_radius * xyz[2];
                normalsArray[ 3 * vertexIndex + 0 ] = xyz[0];
                normalsArray[ 3 * vertexIndex + 1 ] = xyz[1];
                normalsArray[ 3 * vertexIndex + 2 ] = xyz[2];
                uvs_array[ 2 * vertexIndex + 0 ] = u;
                uvs_array[ 2 * vertexIndex + 1 ] = v;
            }
        }
        triangles_array.clear();
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta - 1 ; ++thetaIt ) {
            for( unsigned int phiIt = 0 ; phiIt < nPhi - 1 ; ++phiIt ) {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt+1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt+1) * nTheta;
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuV );
            }
        }
    }


    // Résoudre l'équation d'après le cours
    RaySphereIntersection intersect(const Ray &ray) const {
        
        //TODO calcul l'intersection rayon sphere
        
        float a;
        float b;
        float c;
        float delta2; // delta au carré
        float root1;
        float root2;
        float root_max;
        RaySphereIntersection raySphereIntersection;
        
        // D'après la formule du cours de l'intersection rayon-sphère page 42, on commence par déterminer les coefficients calculables pour exprimer en fonction de t
        a = Vec3::dot(ray.direction(),ray.direction());
        b = 2 * Vec3::dot(ray.direction(),(ray.origin() - this->m_center));
        c = (ray.origin() - this->m_center).squareNorm() - pow(this->m_radius,2);

        delta2 = b * b - (4 * a * c);

        // Puis on récupère l'intersection seulement si le delta est positif
        if(delta2 > 0) {

            root1 = ( -b + sqrt(delta2) ) / ( 2 * a );
            root2 = ( -b - sqrt(delta2) ) / ( 2 * a );

            if(root1 < root2) {

                raySphereIntersection.t = root1;
                root_max = root2;

            } else {

                raySphereIntersection.t = root2;
                root_max = root1;
            }

            if(raySphereIntersection.t < 0.00001) {

                raySphereIntersection.t = FLT_MAX;
                
                return raySphereIntersection;
            } 
            else {

                raySphereIntersection.intersectionExists = true;

                raySphereIntersection.intersection = ray.origin() + raySphereIntersection.t * ray.direction();
                raySphereIntersection.secondintersection = ray.origin() + root_max * ray.direction();

                raySphereIntersection.normal = raySphereIntersection.intersection - this->m_center;
                raySphereIntersection.normal.normalize();
            }
        }

        return raySphereIntersection;
    }
};

#endif
