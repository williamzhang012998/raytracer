

#ifndef CAMERA_H
#define CAMERA_H

#include "Algebra.h"

class Camera {
 public:
  Camera();
  ~Camera();
  void Reset();
  void Orient(Point& eye, Point& focus, Vector& up);
  void Orient(Point& eye, Vector& look, Vector& up);
  void SetViewAngle (double viewAngle);
  void SetNearPlane (double nearPlane);
  void SetFarPlane (double farPlane);
  void SetScreenSize (int screenWidth, int screenHeight);
  void SetCameraType (int cameraType);

  Matrix GetProjectionMatrix();
  Matrix GetModelViewMatrix();

  void RotateV(double angle);
  void RotateU(double angle);
  void RotateW(double angle);
  void Rotate(Point p, Vector axis, double degree);
  void Translate(const Vector &v);

  Point GetEyePoint();
  Vector GetLookVector();
  Vector GetUpVector();
  double GetViewAngle();
  double GetNearPlane();
  double GetFarPlane();
  int GetScreenWidth();
  int GetScreenHeight();

  double GetFilmPlaneDepth();
  double GetScreenWidthRatio();

  Matrix GetRxyz2uvw();
  Matrix GetInvScaleMatrix();


 private:
  void CalculateProjection();
  void CalculateModelview();

  Point m_eye;
  Vector m_look;
  Vector m_up;
  Vector m_u;
  Vector m_v;
  Vector m_w;

  double m_viewAngle;
  double m_clipNear;
  double m_clipFar;
  double m_screenWidth;
  double m_screenHeight;
  int m_cameraType;

  Matrix projection;
  Matrix modelview;

};
#endif

