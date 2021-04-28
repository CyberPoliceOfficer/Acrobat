#ifndef MY_HEADER_H
#define MY_HEADER_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <iostream>
#include <cmath>

using namespace Eigen;

struct SixRSS
{
  float r_b; //Radious of the base
  float r_m; //Radious of the platform
        
  float d_b; //Lenght between base's anchors
  float d_m; //Lenght between platform's anchors
        
  float d; //Platform's arm lenght
  float h; //Servo's arm lenght
        
  float phi; //Angle between servo's arm and platform's base
  float beta;

  VectorXf theta_b;
  VectorXf theta_m;  

  MatrixXf b_k;
  MatrixXf m_k;

  VectorXf beta_k;
  VectorXf phi_k;

  VectorXf u_k;
};

void kinematic_inverse_jacobian (struct SixRSS * , const Ref<const Vector3f> & , const Ref<const Matrix3f> &, Ref<MatrixXf> & );
bool check_pose (struct SixRSS *, const Ref<const Vector3f> &, const Ref<const Matrix3f> &);
void generate_manipulator (struct SixRSS *);
void evaluate_workspace_discretization (struct SixRSS * , float *, float * );
void evaluate_workspace (int , float * , float * , float * );

#endif // MY_HEADER_H