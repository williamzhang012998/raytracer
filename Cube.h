#ifndef CUBE_H
#define CUBE_H

#include "Shape.h"
#include "FlatTNode.h"

class Cube : public Shape {
    public:
        // points = new vertex[vertexCount];
        Point** points;
        // normals
        Vector** normals;
        // triangles = nx3 array of pointers to Points
        Face** triangles;

        Cube() {
            points = NULL;
            normals = NULL;
            triangles = NULL;
        };
        ~Cube() {
        };

        double calculateIntersect(Point eyeObj, Vector rayVObj, Point planeCenter) {
            Vector normal(planeCenter[0], planeCenter[1], planeCenter[2]);
            normal.normalize();
            double t;
            // t = (planeCenter[0] - eyeObj[0]) / rayVObj[0]
            Vector temp_ = planeCenter - eyeObj;
            /* if (!IN_RANGE(dot(rayVObj, normal), 0) || !IN_RANGE(dot(temp_, normal), 0)) { */
            /*   t = (double)dot(temp_, normal) / (double)dot(rayVObj, normal); */
            /* } */
            /* else { */
            /*   t = -1.; */
            /* } */

            if (planeCenter[0] == .5 || planeCenter[0] == -.5) {
                t = (planeCenter[0] - eyeObj[0]) / rayVObj[0];
            }
            else if (planeCenter[1] == .5 || planeCenter[1] == -.5) {
                t = (planeCenter[1] - eyeObj[1]) / rayVObj[1];
            }
            else if (planeCenter[2] == .5 || planeCenter[2] == -.5) {
                t = (planeCenter[2] - eyeObj[2]) / rayVObj[2];
            }
            //double test1 = dot (temp_, normal);
            //Vector _temp_ = normalize(temp_);
            //double test2 = dot (_temp_, normal);
            //if (!(IN_RANGE(test1, test2))) {
            //std::cerr << "test1: " << test1 << " test2: " << test2 << std::endl;
            //}


            Point intersection = eyeObj + rayVObj * t;
            double x = intersection[0];
            double y = intersection[1];
            double z = intersection[2];
            if ((x > -.5 || IN_RANGE(x, -.5)) && (x < .5 || IN_RANGE(x, .5)) && 
                    (y > -.5 || IN_RANGE(y, -.5)) && (y < .5 || IN_RANGE(y, .5)) &&
                    (z > -.5 || IN_RANGE(z, -.5)) && (z < .5 || IN_RANGE(z, .5))) {
                //if ((x > -.5 && x < .5) && 
                // (y > -.5 && y < .5) &&
                //  (z > -.5 && z < .5)) {
                return t;
            }
            else {
                return -1.;
            }
            }

            double Intersect(Point eyePointP, Vector rayV, Matrix transformMatrix) {
                Matrix inv_mat = invert(transformMatrix);
                Point eyeObj = inv_mat * eyePointP;        
                Vector rayVObj = inv_mat * rayV;

                double ts[6];
                ts[0] = calculateIntersect(eyeObj, rayVObj, Point(0.5, 0, 0));
                ts[1] = calculateIntersect(eyeObj, rayVObj, Point(-0.5, 0, 0));
                ts[2] = calculateIntersect(eyeObj, rayVObj, Point(0, 0.5, 0));
                ts[3] = calculateIntersect(eyeObj, rayVObj, Point(0, -0.5, 0));
                ts[4] = calculateIntersect(eyeObj, rayVObj, Point(0, 0, 0.5));
                ts[5] = calculateIntersect(eyeObj, rayVObj, Point(0, 0, -0.5));

                double min = ts[0];
                for (int i = 1; i < 6; i++) {
                    if (ts[i] > 0 || IN_RANGE(ts[i], 0)) {
                        if (IN_RANGE(min, -1.) || (!IN_RANGE(ts[i], -1.) && ts[i] < min)) {
                            min = ts[i];
                        }
                    } 
                }
                //            std::cerr << min << std::endl;
                return min;
            }


            Vector findIsectNormal(Point eyePoint, Vector ray, double dist) {
                Point interPoint(eyePoint + dist * ray);
                //            printf("x: %f y: %f z: %f\n", interPoint.at(0), interPoint.at(1), interPoint.at(2));
                if (IN_RANGE(interPoint[0], .5)) {
                    //              printf("here\n");
                    return Vector(1, 0, 0);
                }
                if (IN_RANGE(interPoint.at(0), -.5)) {
                    return Vector(-1, 0, 0);
                }
                if (IN_RANGE(interPoint.at(1), .5)) {
                    return Vector(0, 1, 0);
                }
                if (IN_RANGE(interPoint.at(1), -.5)) {
                    return Vector(0, -1, 0);
                }
                if (IN_RANGE(interPoint.at(2), .5)) {
                    return Vector(0, 0, 1);
                }
                if (IN_RANGE(interPoint.at(2), -.5)) {
                    return Vector(0, 0, -1);
                }
            }

            void GetTextureColor(Point p, Vector ray, double dist, 
                    float *textureR, float *textureG, float *textureB,
                    FlatTNode *node)
            {
                Point interPoint(p + dist * ray);
                //interPoint = -interPoint;
                int s, t;

                if (IN_RANGE(interPoint.at(0), .5)) {
                    s = (int)(node->textureWidth * node->repeatU * (-interPoint[2] + 0.5)) % node->textureWidth;
                    t = (int)(node->textureHeight * node->repeatV * (-interPoint[1] + 0.5)) % node->textureHeight;
                }
                else if (IN_RANGE(interPoint.at(0), -.5)) {
                    s = (int)(node->textureWidth * node->repeatU * (interPoint[2] + 0.5)) % node->textureWidth;
                    t = (int)(node->textureHeight * node->repeatV * (-interPoint[1] + 0.5)) % node->textureHeight;
                }
                else if (IN_RANGE(interPoint.at(1), .5)) {
                    s = (int)(node->textureWidth * node->repeatU * (interPoint[0] + 0.5)) % node->textureWidth;
                    t = (int)(node->textureHeight * node->repeatV * (interPoint[2] + 0.5)) % node->textureHeight;
                }
                else if (IN_RANGE(interPoint.at(1), -.5)) {
                    s = (int)(node->textureWidth * node->repeatU * (interPoint[0] + 0.5)) % node->textureWidth;
                    t = (int)(node->textureHeight * node->repeatV * (interPoint[2] + 0.5)) % node->textureHeight;
                }
                else if (IN_RANGE(interPoint.at(2), .5)) {
                    s = (int)(node->textureWidth * node->repeatU * (interPoint[0] + 0.5)) % node->textureWidth;
                    t = (int)(node->textureHeight * node->repeatV * (-interPoint[1] + 0.5)) % node->textureHeight;
                }
                else if (IN_RANGE(interPoint.at(2), -.5)) {
                    s = (int)(node->textureWidth * node->repeatU * (-interPoint[0] + 0.5)) % node->textureWidth;
                    t = (int)(node->textureHeight * node->repeatV * (-interPoint[1] + 0.5)) % node->textureHeight;
                } else {
                    cerr << "invalid point" << endl; 
                    cerr << node->GetPrimitive().type << endl;
                    return;
                }

                *textureR = node->texture[t][s].r;
                *textureG = node->texture[t][s].g;
                *textureB = node->texture[t][s].b;
            }

            void draw() {
                int i, j;
                // num triangles on a side
                int numFaces = m_segmentsX * m_segmentsY * 2;
                // int numFaces = m_segmentsX * 2;
                calculatePoints();
                for (j = 0; j < 6; j++) {
                    for (i = 0; i < numFaces; i++) {
                        glBegin(GL_TRIANGLES);
                        normalizeNormal(normals[j][triangles[j][i].p[0]]);
                        glVertex3f(points[j][triangles[j][i].p[0]][0], points[j][triangles[j][i].p[0]][1], points[j][triangles[j][i].p[0]][2]);
                        normalizeNormal(normals[j][triangles[j][i].p[1]]);
                        glVertex3f(points[j][triangles[j][i].p[1]][0], points[j][triangles[j][i].p[1]][1], points[j][triangles[j][i].p[1]][2]);
                        normalizeNormal(normals[j][triangles[j][i].p[2]]);
                        glVertex3f(points[j][triangles[j][i].p[2]][0], points[j][triangles[j][i].p[2]][1], points[j][triangles[j][i].p[2]][2]);
                        glEnd();
                    }
                }
            };


            void drawNormal() {
                int i, j;
                int pointCount = (m_segmentsX + 1) * (m_segmentsY + 1);
                calculatePoints();
                for (j = 0; j < 6; j++) {
                    for (i = 0; i < pointCount; i++) {
                        glBegin(GL_LINES);
                        glVertex3f(points[j][i][0], points[j][i][1], points[j][i][2]);
                        Point endpoint = (normals[j][i] * .1) + points[j][i];
                        glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
                        glEnd();
                    }
                }
            };


            private:
            void calculatePoints () {
                if(m_redraw)
                {
                    int i, j, k;
                    // destroy old arrays
                    if (points) {
                        for (i = 0; i < 6; i++) {
                            delete [] points[i];
                            delete [] triangles[i];
                            delete [] normals[i];
                        }
                        delete [] points;
                        delete [] triangles;
                        delete [] normals;
                    }


                    float lenx, leny;
                    float initX = -.5, initY = .5;

                    // create new arrays
                    int pointCount = (m_segmentsX + 1) * (m_segmentsY + 1);

                    // change 6 to msegmentsX or Y depending on shape
                    points = new Point*[6];
                    triangles = new Face*[6];
                    normals = new Vector*[6];

                    for (i = 0; i < 6; i++) {
                        points[i] = new Point[pointCount];
                        triangles[i] = new Face[m_segmentsX * m_segmentsY * 2];
                        normals[i] = new Vector[pointCount];
                    }

                    lenx = 1. / m_segmentsX;
                    leny = 1. / m_segmentsY;
                    // initializing points
                    // need to do top and bottom faces
                    //        Vector transVector;
                    for (k = 0; k < 6; k++) {
                        for (j =  0; j <= m_segmentsY; j++) {
                            for (i = 0; i <= m_segmentsX; i++) {
                                points[k][(m_segmentsX + 1) * j + i][0] = i * lenx + initX;
                                points[k][(m_segmentsX + 1) * j + i][1] = initY - (j * leny);
                                points[k][(m_segmentsX + 1) * j + i][2] = 0;

                                normals[k][(m_segmentsX + 1) * j + i][0] = 0;
                                normals[k][(m_segmentsX + 1) * j + i][1] = 0;
                                normals[k][(m_segmentsX + 1) * j + i][2] = 1;

                                if (k == 4) {
                                    //                Vector transVector(0, .5, 0);
                                    points[k][(m_segmentsX + 1) * j + i] = rotX_mat(3 * PI / 2) * points[k][(m_segmentsX + 1) * j + i];
                                    points[k][(m_segmentsX + 1) * j + i] = trans_mat(Vector(0, .5, 0)) * points[k][(m_segmentsX + 1) * j + i];

                                    normals[k][(m_segmentsX + 1) * j + i] = rotX_mat(3 * PI / 2) * normals[k][(m_segmentsX + 1) * j + i];
                                }
                                else if (k == 5) {
                                    //Vector transVector(0, -.5, 0);
                                    points[k][(m_segmentsX + 1) * j + i] = rotX_mat(PI / 2) * points[k][(m_segmentsX + 1) * j + i];
                                    points[k][(m_segmentsX + 1) * j + i] = trans_mat(Vector(0, -.5, 0)) * points[k][(m_segmentsX + 1) * j + i];

                                    normals[k][(m_segmentsX + 1) * j + i] = rotX_mat(PI / 2) * normals[k][(m_segmentsX + 1) * j + i];
                                }
                                else {
                                    //Vector transVector(.5 * cos((PI * k / 2) - (PI / 2)), 0, .5 * sin((PI * k / 2) + (PI / 2)));
                                    points[k][(m_segmentsX + 1) * j + i] = rotY_mat(PI * k /2) * points[k][(m_segmentsX + 1) * j + i];
                                    points[k][(m_segmentsX + 1) * j + i] = trans_mat(Vector(.5 * cos((PI * k / 2) - (PI / 2)), 0, .5 * sin((PI * k / 2) + (PI / 2)))) * points[k][(m_segmentsX + 1) * j + i];

                                    normals[k][(m_segmentsX + 1) * j + i] = rotY_mat(PI * k /2) * normals[k][(m_segmentsX + 1) * j + i];
                                }
                            }
                        }
                    }

                    // initializing faces
                    for (k = 0; k < 6; k++) {
                        for (j =  0; j < m_segmentsY; j++) {
                            for (i = 0; i <= m_segmentsX; i++) {
                                // only face to right of point
                                if (i == 0) {
                                    triangles[k][m_segmentsX * j * 2 + i * 2].p[0] = (m_segmentsX + 1) * j + i;
                                    triangles[k][m_segmentsX * j * 2 + i * 2].p[1] = (m_segmentsX + 1) * (j + 1) + i;
                                    triangles[k][m_segmentsX * j * 2 + i * 2].p[2] = (m_segmentsX + 1) * j + (i + 1);
                                }
                                // only face to left of point
                                else if (i == m_segmentsX) {
                                    triangles[k][m_segmentsX * j * 2 + i * 2 - 1].p[0] = (m_segmentsX + 1) * j + i;
                                    triangles[k][m_segmentsX * j * 2 + i * 2 - 1].p[1] = (m_segmentsX + 1) * (j + 1) + (i - 1);
                                    triangles[k][m_segmentsX * j * 2 + i * 2 - 1].p[2] = (m_segmentsX + 1) * (j + 1) + i;   
                                }
                                // face to left then face to right
                                else {
                                    triangles[k][m_segmentsX * j * 2 + i * 2 - 1].p[0] = (m_segmentsX + 1) * j + i;
                                    triangles[k][m_segmentsX * j * 2 + i * 2 - 1].p[1] = (m_segmentsX + 1) * (j + 1) + (i - 1);
                                    triangles[k][m_segmentsX * j * 2 + i * 2 - 1].p[2] = (m_segmentsX + 1) * (j + 1) + i;

                                    triangles[k][m_segmentsX * j * 2 + i * 2].p[0] = (m_segmentsX + 1) * j + i;
                                    triangles[k][m_segmentsX * j * 2 + i * 2].p[1] = (m_segmentsX + 1) * (j + 1) + i;
                                    triangles[k][m_segmentsX * j * 2 + i * 2].p[2] = (m_segmentsX + 1) * j + (i + 1);   
                                }                    
                            }
                        }           
                    }
                    m_redraw = false;
                }    
            };
        };

#endif
