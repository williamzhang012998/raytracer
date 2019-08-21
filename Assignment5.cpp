#define NUM_OPENGL_LIGHTS 8
#define RECUR_THRESHOLD .0001
#include <iostream>
#include <fstream>
#include <string>
#include <GL/glui.h>
#include "Shape.h"
#include "Cube.h"
#include "Cylinder.h"
#include "Cone.h"
#include "Sphere.h"
#include "SceneParser.h"
#include "Camera.h"
#include "FlatTNode.h"
#include <deque>

using namespace std;

/** These are the live variables passed into GLUI ***/
int  isectOnly = 0;
int camRotU = 0;
int camRotV = 0;
int camRotW = 0;
int viewAngle = 45;
float eyeX = 2;
float eyeY = 2;
float eyeZ = 2;
float lookX = -2;
float lookY = -2;
float lookZ = -2;
int max_recursion = 1; 
int antialiasing = 0;

/** These are GLUI control panel objects ***/
int  main_window;
string filenamePath = "data\/tests\/work.xml";
GLUI_EditText* filenameTextField = NULL;
GLubyte* pixels = NULL;
int pixelWidth = 0, pixelHeight = 0;
int screenWidth = 0, screenHeight = 0;

/** these are the global variables used for rendering **/
Cube* cube = new Cube();
Cylinder* cylinder = new Cylinder();
Cone* cone = new Cone();
Sphere* sphere = new Sphere();
SceneParser* parser = NULL;
Camera* camera = new Camera();
std::vector<FlatTNode*> flatTree;
FlatTNode **flatTreeArr;
int flatTreeSize;
void setupCamera();
void updateCamera();
Matrix combineTransformations(
        std::vector<SceneTransformation*> transformations, Matrix m);
void flattenTree(SceneNode* node, std::deque<Matrix> transDeque);
SceneColor rayTrace(Point p, Vector reflected, int iteration, double* intensity, int i, int j, bool& aboveThreshold);
double getMinT(Point p, Vector d, FlatTNode **currMinNode);

void printPV(string label, Point p) {
    cout << label << ": " << p[0] << ", " << p[1] << ", " << p[2] << endl;
}

void printPV(string label, Vector v) {
    cout << label << ": " << v[0] << ", " << v[1] << ", " << v[2] << endl;
}

void setPixel(GLubyte* buf, int x, int y, int r, int g, int b) {
    buf[(y*pixelWidth + x) * 3 + 0] = (GLubyte)r;
    buf[(y*pixelWidth + x) * 3 + 1] = (GLubyte)g;
    buf[(y*pixelWidth + x) * 3 + 2] = (GLubyte)b;
}

void callback_start(int id) {
    cout << "start button clicked!" << endl;

    if (parser == NULL) {
        cout << "no scene loaded yet" << endl;
        return;
    }

    pixelWidth = screenWidth;
    pixelHeight = screenHeight;

    updateCamera();

    if (pixels != NULL) {
        delete pixels;
    }
    pixels = new GLubyte[pixelWidth  * pixelHeight * 3];
    memset(pixels, 0, pixelWidth  * pixelHeight * 3);

    cout << "(w, h): " << pixelWidth << ", " << pixelHeight << endl;

    Point eye_world = camera->GetEyePoint();
    Vector eye_pointV = Vector(eye_world[0], eye_world[1], eye_world[2]);
    eye_pointV.normalize();
    Matrix film2world = inv_trans_mat(-eye_pointV) * transpose(camera->GetRxyz2uvw()) * camera->GetInvScaleMatrix();

    for (int i = 0; i < pixelWidth; i++) {
        for (int j = 0; j < pixelHeight; j++) {
            // Point p_film(((2.0*i) /(double)(pixelWidth - 1)) - 1.0, 
            //         ((2.0*j) /(double)(pixelHeight - 1)) - 1.0,
            //         -1.0);
            Point p_film(((2.0*i) /(double)(pixelWidth)) - 1.0, 
                         ((2.0*j) /(double)(pixelHeight)) - 1.0,
                         -1.0);
            Point p_world = film2world * p_film;
            Vector d(p_world - eye_world);
            d.normalize(); 

            double* intensity = new double[3];
            for (int z = 0; z < 3; z++) {
                intensity[z] = 0.;
            }

            bool aboveThreshold = true;
            SceneColor color1 = rayTrace(eye_world, d, 0, intensity, i, j, aboveThreshold);
            double avgR, avgG, avgB;

            if (antialiasing) {
              p_film = Point(((2.0*i) /(double)(pixelWidth)) - 1.0 - .002, 
                             ((2.0*j) /(double)(pixelHeight)) - 1.0 - .002,
                             -1.0);
              p_world = film2world * p_film;
              d = (p_world - eye_world);
              d.normalize(); 

              for (int z = 0; z < 3; z++) {
                intensity[z] = 0.;
              }
              aboveThreshold = true;
              SceneColor color2 = rayTrace(eye_world, d, 0, intensity, i, j, aboveThreshold);

              p_film = Point(((2.0*i) /(double)(pixelWidth)) - 1.0 + .002, 
                             ((2.0*j) /(double)(pixelHeight)) - 1.0 + .002,
                             -1.0);
              p_world = film2world * p_film;
              d = (p_world - eye_world);
              d.normalize(); 

              for (int z = 0; z < 3; z++) {
                intensity[z] = 0.;
              }
              aboveThreshold = true;
              SceneColor color3 = rayTrace(eye_world, d, 0, intensity, i, j, aboveThreshold);

              p_film = Point(((2.0*i) /(double)(pixelWidth)) - 1.0 - .002, 
                             ((2.0*j) /(double)(pixelHeight)) - 1.0 + .002,
                             -1.0);
              p_world = film2world * p_film;
              d = (p_world - eye_world);
              d.normalize(); 

              for (int z = 0; z < 3; z++) {
                intensity[z] = 0.;
              }
              aboveThreshold = true;
              SceneColor color4 = rayTrace(eye_world, d, 0, intensity, i, j, aboveThreshold);

              p_film = Point(((2.0*i) /(double)(pixelWidth)) - 1.0 + .002, 
                             ((2.0*j) /(double)(pixelHeight)) - 1.0 - .002,
                             -1.0);
              p_world = film2world * p_film;
              d = (p_world - eye_world);
              d.normalize(); 

              for (int z = 0; z < 3; z++) {
                intensity[z] = 0.;
              }
              aboveThreshold = true;
              SceneColor color5 = rayTrace(eye_world, d, 0, intensity, i, j, aboveThreshold);
              avgR = (color1.r + color2.r + color3.r + color4.r + color5.r) / 5.0;
              avgG = (color1.g + color2.g + color3.g + color4.g + color5.g) / 5.0;
              avgB = (color1.b + color2.b + color3.b + color4.b + color5.b) / 5.0;
            }
            else {
              avgR = color1.r;
              avgG = color1.g;
              avgB = color1.b;
            }
                setPixel(pixels, i, j, 
                        avgR * 255., 
                        avgG * 255.,
                        avgB * 255.);

            delete [] intensity;
        }
    }
    std::cout << "Finished!\n";
    glutPostRedisplay();
}

double getMinT(Point p, Vector d, FlatTNode **currMinNode)
{
    double t;
    double min_t;
    for (int k = 0; k < flatTreeSize; k++) {
        switch (flatTreeArr[k]->GetPrimitive().type) {
            case SHAPE_CUBE:
                t = cube->Intersect(p, d, flatTreeArr[k]->GetTMatrix());
                break;
            case SHAPE_SPHERE:
                t = sphere->Intersect(p, d, flatTreeArr[k]->GetTMatrix());
                break;
            case SHAPE_CYLINDER:
                t = cylinder->Intersect(p, d, flatTreeArr[k]->GetTMatrix());
                break;
            case SHAPE_CONE:
                t = cone->Intersect(p, d, flatTreeArr[k]->GetTMatrix());
                break;
            default:
                break;
        }


        if (IN_RANGE(k, 0)) {
            min_t = t;
            *currMinNode = flatTreeArr[k];                    
        } else if (IN_RANGE(min_t, -1) || (!IN_RANGE(t, -1) && t < min_t)) {
            min_t = t;
            *currMinNode = flatTreeArr[k];                    
        }
    }

    return min_t;
}

SceneColor rayTrace(Point p, Vector d, int iteration, double* intensity, int i, int j, bool& aboveThreshold) 
{    
    if (iteration > max_recursion || !aboveThreshold) {
        SceneColor color;
        return color;
    }

    double min_t;
    FlatTNode* currMinNode = NULL;
    min_t = getMinT(p, d, &currMinNode);
 
    if (min_t > 0) {
        if (isectOnly) {
            if (iteration == 0) {
                setPixel(pixels, i, j, 255, 255, 255);
            }
            SceneColor color;
            color.r = 255;
            color.g = 255;
            color.b = 255;
            return color;
        }
        else {
            // convert ray and eye point to object space
            Matrix inv_mat = invert(currMinNode->GetTMatrix());
            Point p_obj(inv_mat * p);
            Vector d_obj(inv_mat * d);

            Vector norm;

            // get the normal in object space
            switch (currMinNode->GetPrimitive().type) {
                case SHAPE_CUBE: 
                    norm = cube->findIsectNormal(p_obj, d_obj, min_t);
                    break;
                case SHAPE_SPHERE:
                    norm = sphere->findIsectNormal(p_obj, d_obj, min_t);
                    break;
                case SHAPE_CYLINDER:
                    norm = cylinder->findIsectNormal(p_obj, d_obj, min_t);
                    break;
                case SHAPE_CONE: 
                    norm = cone->findIsectNormal(p_obj, d_obj, min_t);
                    break;
                default:
                    break;
            }

            // convert normal from object space to world space
            norm = transpose(inv_mat) * norm;
            norm.normalize();

            // get data for illumination calculation
            SceneGlobalData k;
            parser->getGlobalData(k);
            double k_a = k.ka;
            double k_d = k.kd;
            double k_s = k.ks;
            double k_t = k.kt;

            double oal_r = currMinNode->GetPrimitive().material.cAmbient.r;
            double oal_g = currMinNode->GetPrimitive().material.cAmbient.g;
            double oal_b = currMinNode->GetPrimitive().material.cAmbient.b;

            double odl_r = currMinNode->GetPrimitive().material.cDiffuse.r;
            double odl_g = currMinNode->GetPrimitive().material.cDiffuse.g;
            double odl_b = currMinNode->GetPrimitive().material.cDiffuse.b;

            double orl_r = currMinNode->GetPrimitive().material.cReflective.r;
            double orl_g = currMinNode->GetPrimitive().material.cReflective.g;
            double orl_b = currMinNode->GetPrimitive().material.cReflective.b;

            double osl_r = currMinNode->GetPrimitive().material.cSpecular.r;
            double osl_g = currMinNode->GetPrimitive().material.cSpecular.g;
            double osl_b = currMinNode->GetPrimitive().material.cSpecular.b;

            double otl_r = currMinNode->GetPrimitive().material.cTransparent.r;
            double otl_g = currMinNode->GetPrimitive().material.cTransparent.g;
            double otl_b = currMinNode->GetPrimitive().material.cTransparent.b;
            double refractionI = currMinNode->GetPrimitive().material.ior;

            double f = currMinNode->GetPrimitive().material.shininess;

            double intensityR, intensityG, intensityB;

            double sumLightsR = 0;
            double sumLightsG = 0;
            double sumLightsB = 0;

            // get intersection point in world space
            Point intersectionP = p + min_t * d;
            SceneLightData lightData;
            Point epsilonP = intersectionP - 1e-5 * d;
            for (int l = 0; l < parser->getNumLights(); l++) {
                parser->getLightData(l, lightData);
                
                double light_t;
                Vector lightV, lightIntersectV;
                if (lightData.type == LIGHT_POINT) {
                    lightV = lightData.pos - intersectionP;
                    lightIntersectV = lightData.pos - epsilonP;
                    light_t = length(lightIntersectV);
                } else if (lightData.type == LIGHT_DIRECTIONAL) {
                    lightV = -lightData.dir;
                    lightIntersectV = lightV;
                    light_t = length(lightIntersectV);
                } else {
                    cout << "unknown light type" << endl;
                    continue;
                }
                lightV.normalize();
                lightIntersectV.normalize();

                FlatTNode *closestObjNode; 
                double closest_t = getMinT(epsilonP, lightIntersectV, &closestObjNode);

                if (closest_t < 0 || (closest_t > light_t && lightData.type == LIGHT_POINT)) {
                    Vector reflected_light = lightV - 2. * dot(lightV, norm) * norm;
                    reflected_light.normalize();

                    double dot_norm_light = dot(norm, lightV);
                    if (dot_norm_light < 0) {
                        dot_norm_light = 0;
                    }

                    //Vector v(camera->GetLookVector());
                    //v.normalize();
                    //double dot_r_v = dot(reflected_light, v);
                    double dot_r_v = dot(reflected_light, d);
                    if (dot_r_v < 0) {
                        dot_r_v = 0;
                    }

                    sumLightsR += lightData.color.r * (k_d * odl_r * dot_norm_light + k_s * osl_r * pow(dot_r_v, f));
                    sumLightsG += lightData.color.g * (k_d * odl_g * dot_norm_light + k_s * osl_g * pow(dot_r_v, f));
                    sumLightsB += lightData.color.b * (k_d * odl_b * dot_norm_light + k_s * osl_b * pow(dot_r_v, f));
                } else {
                    continue;
                }
            }

            // reflection
            Vector new_d = Vector(d - 2. * dot(d, norm) * norm);
            new_d.normalize();

            rayTrace(epsilonP, new_d, iteration + 1, intensity, i, j, aboveThreshold);

            // refraction
            double refractFrac = 1./refractionI;
            double c2 = sqrt(1. - refractFrac*refractFrac * (1. - SQR(dot(d, norm))));
            if (refractionI > 0 && otl_r > 0 && otl_g > 0 && otl_b > 0) {
              Vector refracted_d = Vector((refractFrac * d) + (refractFrac * -dot(d, norm) - c2) * norm);
              rayTrace(epsilonP, refracted_d, iteration + 1, intensity, i, j, aboveThreshold);
            }

            double recurIntensityR, recurIntensityG, recurIntensityB;
            if (refractionI > 0) {
              recurIntensityR = k_s * orl_r * k_t * otl_r * intensity[0];
              recurIntensityG = k_s * orl_g * k_t * otl_g * intensity[1];
              recurIntensityB = k_s * orl_b * k_t * otl_b * intensity[2];
            }
            else {
              recurIntensityR = k_s * orl_r * intensity[0];
              recurIntensityG = k_s * orl_g * intensity[1];
              recurIntensityB = k_s * orl_b * intensity[2];
            }

            intensityR = k_a * oal_r + sumLightsR + recurIntensityR;
            intensityG = k_a * oal_g + sumLightsG + recurIntensityG;
            intensityB = k_a * oal_b + sumLightsB + recurIntensityB;

            if (currMinNode->hasTexture) {
                float textureR;
                float textureG;
                float textureB;
                float blend = currMinNode->GetPrimitive().material.blend;

                switch (currMinNode->GetPrimitive().type) {
                    case SHAPE_CUBE: 
                        cube->GetTextureColor(p_obj, d_obj, min_t, &textureR, &textureG, &textureB, currMinNode);
                        break;                        
                    case SHAPE_SPHERE:
                      sphere->GetTextureColor(p_obj, d_obj, min_t, &textureR, &textureG, &textureB, currMinNode);
                        break;                        
                    case SHAPE_CYLINDER:
                        cylinder->GetTextureColor(p_obj, d_obj, min_t, &textureR, &textureG, &textureB, currMinNode);
                        break;
                    case SHAPE_CONE: 
                        cone->GetTextureColor(p_obj, d_obj, min_t, &textureR, &textureG, &textureB, currMinNode);
                        break;
                    default:
                        break;
                }

                intensityR = intensityR * (1 - blend) + textureR * blend;
                intensityG = intensityG * (1 - blend) + textureG * blend;
                intensityB = intensityB * (1 - blend) + textureB * blend;
            }

            intensityR = intensityR > 1 ? 1 : intensityR;
            intensity[0] = intensityR;

            intensityG = intensityG > 1 ? 1 : intensityG;
            intensity[1] = intensityG;

            intensityB = intensityB > 1 ? 1 : intensityB;
            intensity[2] = intensityB;

            SceneColor color;
            color.r = intensityR;
            color.g = intensityG;
            color.b = intensityB;

            
            if (intensityR < RECUR_THRESHOLD || intensityG < RECUR_THRESHOLD || intensityB < RECUR_THRESHOLD) {
              aboveThreshold = false;
            }
            return color;

            /*
            if (iteration == 0) {
                setPixel(pixels, i, j, 
                        intensityR * 255., 
                        intensityG * 255.,
                        intensityB * 255.);
            }*/
        }
    }
}



    Matrix combineTransformations(
            std::vector<SceneTransformation*> transformations, Matrix m) {
        int i;
        for (i = 0; i < transformations.size(); i++) {
            if (transformations[i]->type == TRANSFORMATION_TRANSLATE) {
                m = m * trans_mat(transformations[i]->translate);
            }
            else if (transformations[i]->type == TRANSFORMATION_SCALE) {
                m = m * scale_mat(transformations[i]->scale);
            }
            else if (transformations[i]->type == TRANSFORMATION_ROTATE) {
                m = m * rot_mat(transformations[i]->rotate, transformations[i]->angle);
            }
            else {
                m = m * transformations[i]->matrix;
            }
        }
        return m;
    }

    void flattenTree(SceneNode* node, Matrix m) {
        int i;
        int j;

        if(!(node->transformations.empty())) {
            m = combineTransformations(node->transformations, m);
        }

        if (!(node->primitives.empty())) {
            for (i = 0; i < node->primitives.size(); i++) {
                FlatTNode* flatNode = new FlatTNode(*node->primitives[i], m);
                flatTree.push_back(flatNode);
            }
        }

        if (!(node->children.empty())) {
            for (i = 0; i < node->children.size(); i++) {
                flattenTree(node->children[i], m);
            }
        }

        return;
    }

    void callback_load(int id) {
        char curDirName [2048];
        if (filenameTextField == NULL) {
            return;
        }
        printf ("%s\n", filenameTextField->get_text());

        if (parser != NULL) {
            delete parser;
        }
        parser = new SceneParser (filenamePath);
        bool success = parser->parse();
        cout << "success? " << success << endl;
        if (success == false) {
            delete parser;
            vector<FlatTNode*>().swap(flatTree);
            parser = NULL;
        }
        else {
            setupCamera();
            SceneNode* node = parser->getRootNode();
            flatTree.clear();
            Matrix m = Matrix();
            flattenTree(node, m);
            flatTreeArr = &flatTree[0]; 
            flatTreeSize = flatTree.size();
        }
    }


    /***************************************** myGlutIdle() ***********/

    void myGlutIdle(void)
    {
        /* According to the GLUT specification, the current window is
           undefined during an idle callback.  So we need to explicitly change
           it if necessary */
        if (glutGetWindow() != main_window)
            glutSetWindow(main_window);

        glutPostRedisplay();
    }


    /**************************************** myGlutReshape() *************/

    void myGlutReshape(int x, int y)
    {
        float xy_aspect;

        xy_aspect = (float)x / (float)y;
        glViewport(0, 0, x, y);
        camera->SetScreenSize(x, y);

        screenWidth = x;
        screenHeight = y;

        glutPostRedisplay();
    }


    /***************************************** setupCamera() *****************/
    void setupCamera()
    {
        SceneCameraData cameraData;
        parser->getCameraData(cameraData);

        camera->Reset();
        camera->SetViewAngle(cameraData.heightAngle);
        if (cameraData.isDir == true) {
            camera->Orient(cameraData.pos, cameraData.look, cameraData.up);
        }
        else {
            camera->Orient(cameraData.pos, cameraData.lookAt, cameraData.up);
        }

        viewAngle = camera->GetViewAngle();
        Point eyeP = camera->GetEyePoint();
        Vector lookV = camera->GetLookVector();
        eyeX = eyeP[0];
        eyeY = eyeP[1];
        eyeZ = eyeP[2];
        lookX = lookV[0];
        lookY = lookV[1];
        lookZ = lookV[2];
        camRotU = 0;
        camRotV = 0;
        camRotW = 0;
        GLUI_Master.sync_live_all();
    }

    void updateCamera()
    {
        camera->Reset();

        Point guiEye (eyeX, eyeY, eyeZ);
        Point guiLook(lookX, lookY, lookZ);
        camera->SetViewAngle(viewAngle);
        Vector up(camera->GetUpVector());
        camera->Orient(guiEye, guiLook, up);
        camera->RotateU(camRotU);
        camera->RotateV(camRotV);
        camera->RotateW(camRotW);
    }

    /***************************************** myGlutDisplay() *****************/

    void myGlutDisplay(void)
    {
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        if (parser == NULL) {
            return;
        }

        if (pixels == NULL) {
            return;
        }

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glDrawPixels(pixelWidth, pixelHeight, GL_RGB, GL_UNSIGNED_BYTE, pixels);
        glutSwapBuffers();
    }

    void onExit()
    {
        delete cube;
        delete cylinder;
        delete cone;
        delete sphere;
        delete camera;
        if (parser != NULL) {
            delete parser;
        }
        if (pixels != NULL) {
            delete pixels;
        }
    }

    /**************************************** main() ********************/

    int main(int argc, char* argv[])
    {
        atexit(onExit);

        /****************************************/
        /*   Initialize GLUT and create window  */
        /****************************************/

        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
        glutInitWindowPosition(50, 50);
        glutInitWindowSize(500, 500);

        main_window = glutCreateWindow("COMP 175 Assignment 5");
        glutDisplayFunc(myGlutDisplay);
        glutReshapeFunc(myGlutReshape);

        /****************************************/
        /*         Here's the GLUI code         */
        /****************************************/

        GLUI* glui = GLUI_Master.create_glui("GLUI");

        filenameTextField = new GLUI_EditText( glui, "Filename:", filenamePath);
        filenameTextField->set_w(300);
        glui->add_button("Load", 0, callback_load);
        glui->add_button("Start!", 0, callback_start);
        glui->add_checkbox("Isect Only", &isectOnly);
        glui->add_checkbox("Antialiasing", &antialiasing);
        GLUI_Spinner *recursionSpinner = glui->add_spinner("Recursion:", GLUI_SPINNER_INT, &max_recursion);
        recursionSpinner->set_int_limits(0, 10);

        GLUI_Panel *camera_panel = glui->add_panel("Camera");
        (new GLUI_Spinner(camera_panel, "RotateV:", &camRotV))
            ->set_int_limits(-179, 179);
        (new GLUI_Spinner(camera_panel, "RotateU:", &camRotU))
            ->set_int_limits(-179, 179);
        (new GLUI_Spinner(camera_panel, "RotateW:", &camRotW))
            ->set_int_limits(-179, 179);
        (new GLUI_Spinner(camera_panel, "Angle:", &viewAngle))
            ->set_int_limits(1, 179);

        glui->add_column_to_panel(camera_panel, true);

        GLUI_Spinner* eyex_widget = glui->add_spinner_to_panel(camera_panel, "EyeX:", GLUI_SPINNER_FLOAT, &eyeX);
        eyex_widget->set_float_limits(-10, 10);
        GLUI_Spinner* eyey_widget = glui->add_spinner_to_panel(camera_panel, "EyeY:", GLUI_SPINNER_FLOAT, &eyeY);
        eyey_widget->set_float_limits(-10, 10);
        GLUI_Spinner* eyez_widget = glui->add_spinner_to_panel(camera_panel, "EyeZ:", GLUI_SPINNER_FLOAT, &eyeZ);
        eyez_widget->set_float_limits(-10, 10);

        GLUI_Spinner* lookx_widget = glui->add_spinner_to_panel(camera_panel, "LookX:", GLUI_SPINNER_FLOAT, &lookX);
        lookx_widget->set_float_limits(-10, 10);
        GLUI_Spinner* looky_widget = glui->add_spinner_to_panel(camera_panel, "LookY:", GLUI_SPINNER_FLOAT, &lookY);
        looky_widget->set_float_limits(-10, 10);
        GLUI_Spinner* lookz_widget = glui->add_spinner_to_panel(camera_panel, "LookZ:", GLUI_SPINNER_FLOAT, &lookZ);
        lookz_widget->set_float_limits(-10, 10);

        glui->add_button("Quit", 0, (GLUI_Update_CB)exit);

        glui->set_main_gfx_window(main_window);

        /* We register the idle callback with GLUI, *not* with GLUT */
        GLUI_Master.set_glutIdleFunc(myGlutIdle);

        glutMainLoop();

        return EXIT_SUCCESS;
    }



