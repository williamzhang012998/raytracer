#ifndef CYLINDER_H
#define CYLINDER_H

#include "Shape.h"

class Cylinder : public Shape {
    public:
        // points as coordinates
        Point** points;
        // normals
        Vector** normals;
        // triangles array of pointers to Points
        Face** triangles;

        Cylinder() {
            points = NULL;
            normals = NULL;
            triangles = NULL;
        };
        ~Cylinder() {};

        double calcPlaneIsect(Point eyeObj, Vector rayVObj, Point planeCenter) {
            Vector normal(planeCenter[0], planeCenter[1], planeCenter[2]);
            normal.normalize();
            double t;
            /* if (!IN_RANGE(dot(rayVObj, normal), 0) || !IN_RANGE(dot((planeCenter - eyeObj), normal), 0)) { */
            /*   t = (double)dot((planeCenter - eyeObj), normal) / (double)dot(rayVObj, normal); */
            /* } */
            /* else { */
            /*   t = -1.; */
            /* } */
            t = (planeCenter[1] - eyeObj[1]) / rayVObj[1];
            Point intersection = eyeObj + rayVObj * t;
            double x = intersection[0];
            double y = intersection[1];
            double z = intersection[2];            
            if (x*x + z*z < 0.25 || IN_RANGE(x*x + z*z, .25)) {
                return t;
            }
            else {
                return -1.;
            }
        }

        double calcBodyIsect(Point eyeObj, Vector rayVObj) {
            double A = rayVObj.at(0) * rayVObj.at(0) + rayVObj.at(2) * rayVObj.at(2);
            double B = 2*eyeObj.at(0) * rayVObj.at(0) + 2*eyeObj.at(2) * rayVObj.at(2);
            double C = eyeObj.at(0) * eyeObj.at(0) + eyeObj.at(2) * eyeObj.at(2) - 0.25;

            double discriminant = B * B - 4.0 * A * C;

            if (discriminant < 0) {
                return -1;
            } else if (IN_RANGE(discriminant, 0)) {
                return -B / (2.0 * A);
            } else {
                double t1 = (-B + sqrt(discriminant)) / (2.0 * A); 
                double t2 = (-B - sqrt(discriminant)) / (2.0 * A); 

                double t;
                double y;

                if ((t1 > 0 || IN_RANGE(t1, 0)) && (t2 > 0 || IN_RANGE(t2, 0))) {
                    t = (t1 < t2) ? t1 : t2;
                    double o_t = (t1 >= t2) ? t1 : t2;
                    y = eyeObj.at(1) + rayVObj.at(1) * t;
                    if ((y > -0.5 || IN_RANGE(y, -.5)) && (y < 0.5 || IN_RANGE(y, .5))) {
                        return t;
                    } else {
                        y = eyeObj.at(1) + rayVObj.at(1) * o_t;
                        if ((y > -0.5 || IN_RANGE(y, -.5)) && (y < 0.5 || IN_RANGE(y, .5))) {
                            return o_t;
                        } else {
                            return -1.;
                        }
                    }
                } else if (t1 > 0 || IN_RANGE(t1, 0)) {
                    t = t1;
                    y = eyeObj.at(1) + rayVObj.at(1) * t;
                } else if (t2 > 0 || IN_RANGE(t2, 0)) {
                    t = t2;
                    y = eyeObj.at(1) + rayVObj.at(1) * t;
                } else {
                    return -1.;
                }

                if ((y > -0.5 || IN_RANGE(y, -.5)) && (y < 0.5 || IN_RANGE(y, .5))) {
                    return t;
                } else {
                    return -1.;
                }
            }
        }

        double Intersect(Point eyePointP, Vector rayV, Matrix transformMatrix) {
            Matrix inv_mat = invert(transformMatrix);
            Point eyeObj = inv_mat * eyePointP;
            Vector rayVObj = inv_mat * rayV;

            double ts[3];
            ts[0] = calcPlaneIsect(eyeObj, rayVObj, Point(0, 0.5, 0));
            ts[1] = calcPlaneIsect(eyeObj, rayVObj, Point(0, -0.5, 0));
            ts[2] = calcBodyIsect(eyeObj, rayVObj);

            double min = ts[0];
            for (int i = 1; i < 3; i++) {
              if (ts[i] > 0 || IN_RANGE(ts[i], 0)) {
                if (IN_RANGE(min, -1) || (!IN_RANGE(ts[i], -1) && ts[i] < min)) {
                        min = ts[i];
                    }
                } 
            }
            return min;

        }

        Vector findIsectNormal(Point eyePoint, Vector ray, double dist) {
            Point interPoint(eyePoint + ray * dist);
            
            if (IN_RANGE(interPoint.at(1), 0.5)) {
                return Vector(0, 1, 0);
            } else if (IN_RANGE(interPoint.at(1), -0.5)) {
                return Vector(0, -1, 0);
            } else {
              double theta = atan2(interPoint[2], interPoint[0]);
                Vector normal;
                /* if ((interPoint[0] < 0 && interPoint[2] > 0 ) || */
                /*     (interPoint[0] < 0 && interPoint[2] < 0)) { */
                /*   normal = Vector(.5*cos(theta + PI), 0, .5*sin(theta + PI)); */
                /* } */
                /* else if (interPoint[0] > 0 && interPoint[2] < 0) { */
                /*   normal = Vector(.5*cos(theta + 2*PI), 0, .5*sin(theta + 2*PI)); */
                /* } */
                /* else { */
                //while (theta < 0) {
                //  theta += 2. * PI;
                //}
                normal = Vector(cos(theta), 0, sin(theta));
                /* } */
                normal.normalize();
                return normal;
            }
        }

        void GetTextureColor(Point p, Vector ray, double dist, 
                float *textureR, float *textureG, float *textureB,
                FlatTNode *node) 
        {
                Point interPoint(p + dist * ray);
                interPoint = -interPoint;
                int s, t;

                if (IN_RANGE(interPoint.at(1), .5)) {
                    s = (int)(node->textureWidth * node->repeatU * (interPoint[0] + 0.5)) % node->textureWidth;
                    t = (int)(node->textureHeight * node->repeatV * (-interPoint[2] + 0.5)) % node->textureHeight;
                }
                else if (IN_RANGE(interPoint.at(1), -.5)) {
                    s = (int)(node->textureWidth * node->repeatU * (interPoint[0] + 0.5)) % node->textureWidth;
                    t = (int)(node->textureHeight * node->repeatV * (interPoint[2] + 0.5)) % node->textureHeight;
                } else {
                    double theta = -atan2(interPoint[2], interPoint[0]);

                    while (theta < 0) {
                        theta += 2*PI;
                    }

                    s = (int)(node->textureWidth * node->repeatU * (theta / (2 * PI))) % node->textureWidth;
                    t = (int)(node->textureHeight * node->repeatV * (interPoint[1] + 0.5)) % node->textureHeight;
                }

                *textureR = node->texture[t][s].r;
                *textureG = node->texture[t][s].g;
                *textureB = node->texture[t][s].b;
        }

        void draw() {
            calculatePoints();

            int i, j;
            int sideCount = m_segmentsX;
            int numFaces = 2 * m_segmentsY;
            for (j = 0; j < sideCount; j++) {
                // draw side triangles
                for (i = 0; i < numFaces; i++) {
                    glBegin(GL_TRIANGLES);
                    normalizeNormal(normals[j][triangles[j][i].p[0]]);
                    glVertex3f(
                            points[j][triangles[j][i].p[0]][0], 
                            points[j][triangles[j][i].p[0]][1], 
                            points[j][triangles[j][i].p[0]][2]);
                    normalizeNormal(normals[j][triangles[j][i].p[1]]);
                    glVertex3f(
                            points[j][triangles[j][i].p[1]][0], 
                            points[j][triangles[j][i].p[1]][1], 
                            points[j][triangles[j][i].p[1]][2]);
                    normalizeNormal(normals[j][triangles[j][i].p[2]]);
                    glVertex3f(
                            points[j][triangles[j][i].p[2]][0], 
                            points[j][triangles[j][i].p[2]][1], 
                            points[j][triangles[j][i].p[2]][2]);
                    glEnd();
                }
                // draw top triangles
                glBegin(GL_TRIANGLES);
                normalizeNormal(
                        normals[sideCount][triangles[sideCount][j].p[0]]);
                glVertex3f(
                        points[sideCount][triangles[sideCount][j].p[0]][0],
                        points[sideCount][triangles[sideCount][j].p[0]][1],
                        points[sideCount][triangles[sideCount][j].p[0]][2]);
                normalizeNormal(
                        normals[sideCount][triangles[sideCount][j].p[1]]);
                glVertex3f(
                        points[sideCount][triangles[sideCount][j].p[1]][0],
                        points[sideCount][triangles[sideCount][j].p[1]][1],
                        points[sideCount][triangles[sideCount][j].p[1]][2]);
                normalizeNormal(
                        normals[sideCount][triangles[sideCount][j].p[2]]);
                glVertex3f(
                        points[sideCount][triangles[sideCount][j].p[2]][0], 
                        points[sideCount][triangles[sideCount][j].p[2]][1], 
                        points[sideCount][triangles[sideCount][j].p[2]][2]);
                glEnd();
                // draw bottom triangles
                glBegin(GL_TRIANGLES);
                normalizeNormal(
                        normals[sideCount + 1][triangles[sideCount][j].p[0]]);
                glVertex3f(
                        points[sideCount + 1][triangles[sideCount][j].p[0]][0],
                        points[sideCount + 1][triangles[sideCount][j].p[0]][1],
                        points[sideCount + 1][triangles[sideCount][j].p[0]][2]);
                normalizeNormal(
                        normals[sideCount + 1][triangles[sideCount][j].p[1]]);
                glVertex3f(
                        points[sideCount + 1][triangles[sideCount][j].p[1]][0],
                        points[sideCount + 1][triangles[sideCount][j].p[1]][1],
                        points[sideCount + 1][triangles[sideCount][j].p[1]][2]);
                normalizeNormal(
                        normals[sideCount + 1][triangles[sideCount][j].p[2]]);
                glVertex3f(
                        points[sideCount + 1][triangles[sideCount][j].p[2]][0], 
                        points[sideCount + 1][triangles[sideCount][j].p[2]][1], 
                        points[sideCount + 1][triangles[sideCount][j].p[2]][2]);
                glEnd();
            }
        };

        void drawNormal() {
            calculatePoints();

            int i, j;
            int sideCount = m_segmentsX;
            int sidePointCount = 2 * (m_segmentsY + 1);
            Point endpoint;
            for (j = 0; j < sideCount; j++) {
                // draw side normals
                for (i = 0; i < sidePointCount; i++) {
                    glBegin(GL_LINES);
                    glVertex3f(
                            points[j][i][0], 
                            points[j][i][1], 
                            points[j][i][2]);
                    endpoint = (normals[j][i] * .1) + points[j][i];
                    glVertex3f(
                            endpoint[0], 
                            endpoint[1], 
                            endpoint[2]);
                    glEnd();
                }
            }
            for (j = 0; j <= sideCount + 1; j++) {
                // draw top normals
                glBegin(GL_LINES);
                glVertex3f(
                        points[sideCount][j][0], 
                        points[sideCount][j][1], 
                        points[sideCount][j][2]);
                endpoint = 
                    (normals[sideCount][j] * .1) + 
                    points[sideCount][j];
                glVertex3f(
                        endpoint[0], 
                        endpoint[1], 
                        endpoint[2]);
                glEnd();
                // draw bottom normals
                glBegin(GL_LINES);
                glVertex3f(
                        points[sideCount + 1][j][0], 
                        points[sideCount + 1][j][1], 
                        points[sideCount + 1][j][2]);
                endpoint = 
                    (normals[sideCount + 1][j] * .1) + 
                    points[sideCount + 1][j];
                glVertex3f(
                        endpoint[0], 
                        endpoint[1], 
                        endpoint[2]);
                glEnd();
            }
        };

        void calculatePoints() {
            if(m_redraw)
            {
                int i, j, k;
                float initY;
                float lenY;
                float segAngle;

                if (points) {
                    for (i = 0; i < prev_segmentsX + 2; i++) {
                        delete [] points[i];
                        delete [] triangles[i];
                        delete [] normals[i];
                    }
                    delete [] points;
                    delete [] triangles;
                    delete [] normals;
                }

                int sideCount = m_segmentsX;
                int sidePointCount = 2 * (m_segmentsY + 1);

                segAngle = 2 * PI / sideCount;

                // length of triangle segements in a single face
                initY = 0.5;
                lenY = 1. / m_segmentsY;

                // create new arrays
                // additional 2 faces for the top and bottom of the cylinder
                points = new Point*[sideCount + 2];
                normals = new Vector*[sideCount + 2];
                triangles = new Face*[sideCount + 2];

                // init arrays for side points 
                for (i = 0; i < sideCount; i++) {
                    points[i] = new Point[sidePointCount];
                    normals[i] = new Vector[sidePointCount];
                    triangles[i] = new Face[2 * m_segmentsY];
                }
                // init arrays for top/bottom points
                for (i = 0; i < 2; i++) {
                    // additional 2 points for: 
                    // + 1 for ending point that is the same as the start point
                    // + 1 for origin point in center of top and bottom 
                    points[sideCount + i] = new Point[m_segmentsX + 2];
                    normals[sideCount + i] = new Vector[m_segmentsX + 2];
                    triangles[sideCount + i] = new Face[m_segmentsX];
                }

                // initialize points for each side
                Vector transVector;
                for ( k = 0; k < sideCount; k++ ){
                    for ( i = 0; i < 2; i++ ){
                        for ( j = 0; j <= m_segmentsY; j++ ){
                            // side points
                            points[k][( 2 * j ) + i][0] = 
                                0.5 * cos((segAngle * i) + (k * segAngle)); 
                            points[k][( 2 * j ) + i][1] =
                                initY - ( lenY * j ); 
                            points[k][( 2 * j ) + i][2] = 
                                0.5 * sin(segAngle * i + (k * segAngle));  

                            normals[k][( 2 * j ) + i][0] = 
                                cos((segAngle * i) + (k * segAngle)); 
                            normals[k][( 2 * j ) + i][1] = 0;
                            normals[k][( 2 * j ) + i][2] =
                                sin(segAngle * i + (k * segAngle));  
                        }
                    } 
                }
                for ( k = 0; k <= sideCount; k++ ){
                    // top points
                    points[sideCount][k][0] = 
                        0.5 * cos(k * segAngle); 
                    points[sideCount][k][1] = 0;
                    points[sideCount][k][2] = 
                        0.5 * sin(k * segAngle);  
                    transVector = Vector(0, 0.5, 0);
                    points[sideCount][k] = 
                        trans_mat(transVector) * 
                        points[sideCount][k]; 
                    normals[sideCount][k][0] = 0;
                    normals[sideCount][k][1] = 1;
                    normals[sideCount][k][2] = 0;
                    // bottom points
                    points[sideCount + 1][k][0] = 
                        0.5 * cos(k * segAngle); 
                    points[sideCount + 1][k][1] = 0;
                    points[sideCount + 1][k][2] = 
                        0.5 * sin(k * segAngle);  
                    transVector = Vector(0, -0.5, 0);
                    points[sideCount + 1][k] = 
                        trans_mat(transVector) * 
                        points[sideCount + 1][k]; 
                    normals[sideCount + 1][k][0] = 0;
                    normals[sideCount + 1][k][1] = -1;
                    normals[sideCount + 1][k][2] = 0;
                }
                // add origin points for top/bottom
                points[sideCount][m_segmentsX + 1][0] = 0; 
                points[sideCount][m_segmentsX + 1][1] = 0;  
                points[sideCount][m_segmentsX + 1][2] = 0; 
                transVector = Vector(0, 0.5, 0);
                points[sideCount][m_segmentsX + 1] = 
                    trans_mat(transVector) * 
                    points[sideCount][m_segmentsX + 1];
                normals[sideCount][m_segmentsX + 1][0] = 0;
                normals[sideCount][m_segmentsX + 1][1] = 1;
                normals[sideCount][m_segmentsX + 1][2] = 0;

                points[sideCount + 1][m_segmentsX + 1][0] = 0; 
                points[sideCount + 1][m_segmentsX + 1][1] = 0;  
                points[sideCount + 1][m_segmentsX + 1][2] = 0; 
                transVector = Vector(0, -0.5, 0);
                points[sideCount + 1][m_segmentsX + 1] = 
                    trans_mat(transVector) * 
                    points[sideCount + 1][m_segmentsX + 1];
                normals[sideCount + 1][m_segmentsX + 1][0] = 0;
                normals[sideCount + 1][m_segmentsX + 1][1] = -1;
                normals[sideCount + 1][m_segmentsX + 1][2] = 0;

                // group triangles into side faces
                for ( k = 0; k < sideCount; k++ ){
                    for ( j = 0; j < m_segmentsY; j++ ){
                        // left triangle
                        triangles[k][( 2 * j )].p[0] = 
                            ( 2 * j );
                        triangles[k][( 2 * j )].p[1] = 
                            ( 2 * (j + 1));
                        triangles[k][( 2 * j )].p[2] = 
                            ( 2 * (j + 1)) + 1;
                        // right triangle
                        triangles[k][( 2 * j ) + 1].p[0] = 
                            ( 2 * j ) + 1;
                        triangles[k][( 2 * j ) + 1].p[1] = 
                            ( 2 * j );
                        triangles[k][( 2 * j ) + 1].p[2] = 
                            ( 2 * (j + 1)) + 1;
                    }
                }
                // group triangles into top/bottom faces
                for ( k = 0; k < 2; k++ ){
                    for ( i = 0; i < m_segmentsX; i++ ){
                        // origin point
                        triangles[sideCount + k][i].p[0] = m_segmentsX + 1;
                        triangles[sideCount + k][i].p[1] = i;
                        triangles[sideCount + k][i].p[2] = i + 1;
                    }
                }

                m_redraw = false;
            }
        };

};
#endif
