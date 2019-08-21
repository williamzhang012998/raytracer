#include "Camera.h"

Camera::Camera() {
  m_eye = Point();
  m_look = Vector();
  m_up = Vector();
  m_viewAngle = PI / 4.;
  m_clipNear = 0 + EPSILON;
  m_clipFar = 100.;
  m_screenWidth = 0.;
  m_screenHeight = 0.;
  m_cameraType = 1;
  projection = Matrix();
  modelview = Matrix();
}

Camera::~Camera() {
}

void Camera::Reset() {
  projection = Matrix();
  modelview = Matrix();
  m_viewAngle = PI / 4.;
  m_clipNear = 0 + EPSILON;
  m_clipFar = 100.;
  m_cameraType = 1; 
}

void Camera::Orient(Point& eye, Point& focus, Vector& up) {
  m_look = focus - eye;
  Orient(eye, m_look, up);
}

void Camera::Orient(Point& eye, Vector& look, Vector& up) {
  m_eye = eye;
  m_look = look;
  m_up = up;

  m_w = -1.0 * m_look;
  m_w.normalize();
  m_u = cross(m_up, m_w);
  m_u.normalize();
  m_v = cross(m_w, m_u);
  CalculateProjection();
  CalculateModelview();
}

void Camera::CalculateProjection() {
  /***** for perspective camera *****/
  // M_pp * S_uvw
  double aspect_ratio = m_screenHeight/m_screenWidth;
  Matrix s_uvw = Matrix(
                        1./(tan(m_viewAngle/2.)*m_clipFar*aspect_ratio), 0, 0, 0,
                        0, 1./(tan(m_viewAngle/2.)*m_clipFar*aspect_ratio), 0, 0,
                        0, 0, 1./m_clipFar, 0,
                        0, 0, 0, 1);
  double c = -m_clipNear / m_clipFar;
  Matrix m_pp = Matrix(
                       1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, -1./(c+1.), c/(c+1.),
                       0, 0, -1, 0);
  
  /***** for orthographic camera *****/
  Matrix s_orth = Matrix(
                         2./2., 0, 0, 0,
                         0, 2./2., 0, 0, 
                         0, 0, 1./(m_clipNear - m_clipFar), 0,
                         0, 0, 0, 1
                         );


  if (m_cameraType == 1) {
    projection = m_pp * s_uvw;
  }
  else {
    projection = s_orth;
  }  
}

Matrix Camera::GetProjectionMatrix() {   
  return projection;
}

void Camera::SetCameraType (int cameraType) {
  m_cameraType = cameraType;
  CalculateProjection();
  CalculateModelview();
}

void Camera::SetViewAngle (double viewAngle) {
  m_viewAngle = DEG_TO_RAD(viewAngle);
  CalculateProjection();
}

void Camera::SetNearPlane (double nearPlane) {
  m_clipNear = nearPlane;
  CalculateProjection();
}

void Camera::SetFarPlane (double farPlane) {
  m_clipFar = farPlane;
  CalculateProjection();
}

void Camera::SetScreenSize (int screenWidth, int screenHeight) {
  m_screenWidth = screenWidth;
  m_screenHeight = screenHeight;
  CalculateProjection();
}

void Camera::CalculateModelview() {
  /***** for perspective camera *****/
  // R_xyz2uvw * T_xyz
  Matrix t = Matrix(
                    1, 0, 0, -m_eye.at(0),
                    0, 1, 0, -m_eye.at(1),
                    0, 0, 1, -m_eye.at(2),
                    0, 0, 0, 1);
  
  /***** for orthographic camera *****/
  Point parallelEye = Point(m_eye + m_clipNear * m_look);
  Matrix tPara = Matrix(
                        1, 0, 0, -parallelEye.at(0),
                        0, 1, 0, -parallelEye.at(1),
                        0, 0, 1, -parallelEye.at(2),
                        0, 0, 0, 1);            

  Matrix r = Matrix (
                    m_u.at(0), m_u.at(1), m_u.at(2), 0,
                    m_v.at(0), m_v.at(1), m_v.at(2), 0,
                    m_w.at(0), m_w.at(1), m_w.at(2), 0,
                    0, 0, 0, 1);

  // if perspective camera
  if (m_cameraType == 1) {
    modelview = r * t;
  }
  else {
    modelview = r * tPara;
  }
}

Matrix Camera::GetRxyz2uvw() {
  Matrix r = Matrix (
                    m_u.at(0), m_u.at(1), m_u.at(2), 0,
                    m_v.at(0), m_v.at(1), m_v.at(2), 0,
                    m_w.at(0), m_w.at(1), m_w.at(2), 0,
                    0, 0, 0, 1);
  return r;
}

Matrix Camera::GetInvScaleMatrix() {
  double aspect_ratio = m_screenHeight/m_screenWidth;
  // Matrix s_uvw = Matrix(
  //                       (tan(m_viewAngle/2.)*m_clipFar), 0, 0, 0,
  //                       0, (tan(m_viewAngle/2.)*m_clipFar*aspect_ratio), 0, 0,
  //                       0, 0, m_clipFar, 0,
  //                       0, 0, 0, 1);
  Matrix s_uvw = Matrix(
                        1./(1./(tan(m_viewAngle/2.)*m_clipFar*aspect_ratio)), 0, 0, 0,
                        0, 1./(1./(tan(m_viewAngle/2.)*m_clipFar*aspect_ratio)), 0, 0,
                        0, 0, 1./(1./m_clipFar), 0,
                        0, 0, 0, 1);

  return s_uvw;

}

Matrix Camera::GetModelViewMatrix() {
  return modelview;
}

// up-down rotation
void Camera::RotateU(double angle) {
  Matrix mu = rotX_mat(DEG_TO_RAD(angle));
  modelview = mu * modelview;
}

// left-right rotation
void Camera::RotateV(double angle) {
  Matrix mv = rotY_mat(DEG_TO_RAD(angle));
  modelview = mv * modelview;
}

// tilt rotation
void Camera::RotateW(double angle) {
  Matrix mw = rotZ_mat(DEG_TO_RAD(angle));
  modelview = mw * modelview;
}

void Camera::Translate(const Vector &v) {
  // translate eye by v
}

void Camera::Rotate(Point p, Vector axis, double degrees) {
  // ? 
}

Point Camera::GetEyePoint() {
  return m_eye;
}

Vector Camera::GetLookVector() {
  return m_look;
}

Vector Camera::GetUpVector() {
  return m_up;
}

double Camera::GetViewAngle() {
  return RAD_TO_DEG(m_viewAngle);
}

double Camera::GetNearPlane() {
  return m_clipNear;
}

double Camera::GetFarPlane() {
  return m_clipFar;
}

int Camera::GetScreenWidth() {
  return m_screenWidth;
}

int Camera::GetScreenHeight() {
  return m_screenHeight;
}

double Camera::GetFilmPlaneDepth() {
  return m_clipFar - m_clipNear;
}

double Camera::GetScreenWidthRatio() {
  return m_screenWidth/m_screenHeight;
}
