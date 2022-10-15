#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"

#include <GL/glut.h>

#define MESH_INTERSECTION 0
#define SPHERE_INTERSECTION 1
#define SQUARE_INTERSECTION 2

enum LightType {
    LightType_Spherical,
    LightType_Quad
};


struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Mesh quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

};

struct RaySceneIntersection{
    bool intersectionExists;
    unsigned int typeOfIntersectedObject;
    unsigned int objectIndex;
    float t;
    RayTriangleIntersection rayMeshIntersection; // Index 0
    RaySphereIntersection raySphereIntersection; // Index 1
    RaySquareIntersection raySquareIntersection; // Index 2
    RaySceneIntersection() : intersectionExists(false) , t(FLT_MAX) {} // Si à la fin du calcul t est toujours égal à FLT_MAX alors il n'y a pas eu d'intersection
};



class Scene {
    std::vector< Mesh > meshes;
    std::vector< Sphere > spheres;
    std::vector< Square > squares;
    std::vector< Light > lights;

public:


    Scene() {
    }

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for( unsigned int It = 0 ; It < meshes.size() ; ++It ) {
            Mesh const & mesh = meshes[It];
            mesh.draw();
        }
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            sphere.draw();
        }
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) {
            Square const & square = squares[It];
            square.draw();
        }
    }


    RaySceneIntersection computeIntersection(Ray const & ray) {
        
        //TODO calculer les intersections avec les objets de la scene et garder la plus proche
        
        int nb_spheres;
        int nb_squares;
        int nb_meshes;
        RaySceneIntersection result;

        nb_meshes = this->meshes.size();

        // On détermine l'intersection du rayon pour chaque mesh et on sauvegarde la plus proche

        for(int i = 0; i < nb_meshes; i++) {

            RayTriangleIntersection rayMeshIntersection = this->meshes[i].intersect(ray);

            if(rayMeshIntersection.intersectionExists) {

                if(rayMeshIntersection.t < result.t) {

                    result.intersectionExists = true;
                    result.objectIndex = i;
                    result.typeOfIntersectedObject = MESH_INTERSECTION;
                    result.rayMeshIntersection = rayMeshIntersection;
                    result.t = rayMeshIntersection.t;
                }
            }
        }

        nb_spheres = this->spheres.size();

        // On détermine l'intersection du rayon pour chaque sphère et on sauvegarde la plus proche

        for(int i = 0; i < nb_spheres; i++) {

            RaySphereIntersection raySphereIntersection = this->spheres[i].intersect(ray);

            if(raySphereIntersection.intersectionExists) {

                if(raySphereIntersection.t < result.t) {

                    result.intersectionExists = true;
                    result.objectIndex = i;
                    result.typeOfIntersectedObject = SPHERE_INTERSECTION;
                    result.raySphereIntersection = raySphereIntersection;
                    result.t = raySphereIntersection.t;
                }
            }
        }

        nb_squares = this->squares.size();

        // On détermine l'intersection du rayon pour chaque quad et on sauvegarde la plus proche

        for(int i = 0; i < nb_squares; i++) {

            RaySquareIntersection raySquareIntersection = this->squares[i].intersect(ray);

            if(raySquareIntersection.intersectionExists) {

                if(raySquareIntersection.t < result.t) {

                    result.intersectionExists = true;
                    result.objectIndex = i;
                    result.typeOfIntersectedObject = SQUARE_INTERSECTION;
                    result.raySquareIntersection = raySquareIntersection;
                    result.t = raySquareIntersection.t;
                }
            }
        }

        return result;
    }

    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces ) {

        // TODO RaySceneIntersection raySceneIntersection = computeIntersection(ray);

        RaySceneIntersection raySceneIntersection;
        Mesh intersectedMesh;
        Vec3 color = Vec3(0.,0.,0.);
        Vec3 currentIntersection;
        RaySceneIntersection shadowIntersection;
        Ray shadowRay;
        float shadowSoftness;

        raySceneIntersection = computeIntersection(ray);
                
        if(raySceneIntersection.intersectionExists) {

            // Modèle de Phong

            Vec3 I; // Intensité de la lumière 
            Vec3 Ka; // Réflexion ambiante de l'objet
            Vec3 Kd; // Réflexion diffuse de l'objet
            Vec3 Ks; // Réflexion spéculaire de l'objet

            Vec3 n; // Normale à la surface de l'objet
            Vec3 l; // Direction de la lumière
            Vec3 r; // Direction de réflexion de la lumière sur un miroir
            Vec3 v; // Direction de la surface de l'objet vers la source du rayon
            float m; // Brillance de l'objet

            float ambient;
            float diffuse;
            float specular;

            switch(raySceneIntersection.typeOfIntersectedObject) {

                case MESH_INTERSECTION:

                    currentIntersection = raySceneIntersection.rayMeshIntersection.intersection;

                    intersectedMesh = this->meshes[raySceneIntersection.objectIndex];

                    n = raySceneIntersection.rayMeshIntersection.normal;
                    n.normalize();

                    break;

                case SPHERE_INTERSECTION:

                    currentIntersection = raySceneIntersection.raySphereIntersection.intersection;

                    intersectedMesh = this->spheres[raySceneIntersection.objectIndex];

                    n = raySceneIntersection.raySphereIntersection.normal;
                    n.normalize();

                    break;

                case SQUARE_INTERSECTION: 

                    currentIntersection = raySceneIntersection.raySquareIntersection.intersection;

                    intersectedMesh = this->squares[raySceneIntersection.objectIndex];

                    n = raySceneIntersection.raySquareIntersection.normal;
                    n.normalize();

                    break;
            }

            unsigned int nb_light = this->lights.size();

            for(unsigned int i = 0; i < nb_light; i++) {

                l = this->lights[i].pos - currentIntersection;
                l.normalize();

                r  = 2 * Vec3::dot(l,n) * n - l;
                r.normalize();

                v = ray.origin() - currentIntersection;
                v.normalize();
            
                switch(intersectedMesh.material.type) {

                    case Material_Diffuse_Blinn_Phong:

                        l = this->lights[i].pos - currentIntersection;
                        l.normalize();

                        r  = 2 * Vec3::dot(l,n) * n - l;
                        r.normalize();

                        v = ray.origin() - currentIntersection;
                        v.normalize();

                        m = intersectedMesh.material.shininess;
                        
                        for(unsigned int j = 0; j < 3; j++) {

                            I[j] = this->lights[i].material[j];
                            Ka[j] = intersectedMesh.material.ambient_material[j];
                            Kd[j] = intersectedMesh.material.diffuse_material[j];
                            Ks[j] = intersectedMesh.material.specular_material[j];

                            ambient = Ka[j] * I[j];
                            diffuse = Kd[j] * I[j] * Vec3::dot(n,l);
                            specular = Ks[j] * I[j] * pow(std::max(Vec3::dot(r,v),(float)0.),m);

                            // Pour chaque composante RGB on calcule son ombrage
                            color[j] = ambient + diffuse + specular;
                            
                            if(color[j] < 0.00001) {

                                color[j] = 0;
                            }
                        }

                        break;

                    case Material_Mirror:

                        if(NRemainingBounces > 0)  {

                            l = ray.direction();
                            l.normalize();

                            r =  l - 2 * Vec3::dot(l,n) * n;
                            r.normalize();

                            Ray newRay = Ray(currentIntersection,r);
                            color = rayTraceRecursive(newRay,NRemainingBounces-1);
                        }

                        break;

                    case Material_Glass:

                        break;
                }
            }

            // Envoi d'un rayon d'ombre depuis le point d'intersection vers chaque lumière (ombre dure) (désactivé)
            // for(unsigned int i = 0; i < nb_light; i++) {

            //     l = this->lights[i].pos - currentIntersection;

            //     shadowRay = Ray(currentIntersection,l);     
                
            //     shadowIntersection = computeIntersection(shadowRay);

            //     if(shadowIntersection.t > 0.00001 && shadowIntersection.t < 1.) {
            //         color = Vec3(0.,0.,0.);
            //     }
            // }

            // Envoi d'un rayon d'ombre depuis le point d'intersection vers chaque lumière (ombre douce)
            int nb_rays = 10;

            for(unsigned int i = 0; i < nb_light; i++) {

                shadowSoftness = computeSoftness(i,nb_rays,currentIntersection);

                color *= 1 - shadowSoftness;
            }

        } 

        return color;
    }

    float computeSoftness(int i, int nb_rays, Vec3 intersection) {

        float horizontal_step;
        float vertical_step;
        float radius;
        Vec3 ray_pos_origin;
        Vec3 light_vector;
        Ray shadow_ray;
        float t;
        
        radius = lights[i].radius / 2;
        ray_pos_origin = Vec3(lights[i].pos[0] - radius, lights[i].pos[1] - radius, lights[i].pos[2] - radius);
        
        unsigned int shadow_count = 0;

        for(int i = 0; i < nb_rays; i++) {

            horizontal_step = (float)((std::rand() / (float)(RAND_MAX / (2 * radius))));
            vertical_step = (float)((std::rand() / (float)(RAND_MAX / (2 * radius))));

            light_vector = Vec3(ray_pos_origin[0] + horizontal_step,ray_pos_origin[1],ray_pos_origin[2] + vertical_step) - intersection;
            light_vector.normalize();

            shadow_ray = Ray(intersection,light_vector);

            t = closestIntersection(shadow_ray);

            if(t > 0.00001 && t < 1) {

                shadow_count++;
            } 

        }

        return (float)shadow_count / nb_rays;
    }

    float closestIntersection(const Ray& ray) {

        RayTriangleIntersection rayTriangleIntersection;
        RaySphereIntersection raySphereIntersection;
        RaySquareIntersection raySquareIntersection;

        unsigned int mesh_count = meshes.size();
        unsigned int sphere_count = spheres.size();
        unsigned int square_count = squares.size();


        for(unsigned int i = 0; i < mesh_count; i++) {

            rayTriangleIntersection = meshes[i].intersect(ray);

            if(rayTriangleIntersection.intersectionExists) {

                return rayTriangleIntersection.t;
            }
        }

        for(unsigned int i = 0; i < sphere_count; i++) {

            raySphereIntersection = spheres[i].intersect(ray);

            if(raySphereIntersection.intersectionExists) {

                return raySphereIntersection.t;
            }
        }

        for(unsigned int i = 0; i < square_count; i++) {

            raySquareIntersection = squares[i].intersect(ray);

            if(raySquareIntersection.intersectionExists) {

                return raySquareIntersection.t;
            }
        }

        return FLT_MAX;
    }

    //TODO appeler la fonction recursive
    Vec3 rayTrace( Ray const & rayStart ) {

        Vec3 color = rayTraceRecursive(rayStart,1);

        return color;
    }

    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1. , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1. , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.,1.,0. );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
    }

    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.8,0.8,0.8 );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
    }

    void setup_cornell_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.25,0.25,0.25 );
            s.material.specular_material = Vec3( 0.25,0.25,0.25 );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.5,0.5,0.5 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.5,0.5,0.5 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        // { //Front Wall
        //     squares.resize( squares.size() + 1 );
        //     Square & s = squares[squares.size() - 1];
        //     s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
        //     s.translate(Vec3(0., 0., -2.));
        //     s.scale(Vec3(2., 2., 1.));
        //     s.rotate_y(180);
        //     s.build_arrays();
        //     s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
        //     s.material.specular_material = Vec3( 1.0,1.0,1.0 );
        //     s.material.shininess = 16;
        // }


        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.,1.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

    }

    void setup_mesh_in_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            meshes.resize(spheres.size() + 1);
            Mesh &m = meshes[meshes.size() - 1];
            m.loadOFF("data/icosahedron.off");
            m.centerAndScaleToUnit();
            m.translate(Vec3(-0.5,-1.,-1.));
            m.build_arrays();
            m.material.type = Material_Mirror;
            m.material.diffuse_material = Vec3(1., 0., 0.);
            m.material.specular_material = Vec3(1., 1., 1.);
            m.material.shininess = 16;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.25,0.25,0.25 );
            s.material.specular_material = Vec3( 0.25,0.25,0.25 );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.5,0.5,0.5 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.5,0.5,0.5 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }
    }
};

#endif
