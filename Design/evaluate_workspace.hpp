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

  float joint_limit;

  Vector<float, 6> theta_b;
  Vector<float, 6> theta_m;  

  Matrix<float, 3, 6> b_k;
  Matrix<float, 3, 6> m_k;

  Vector<float, 6> beta_k;
  Vector<float, 6> phi_k;

  Matrix<float, 3, 6> u_k;
};

void get_h_k (struct SixRSS * , const Ref<const VectorXf>& , Ref<MatrixXf> );
void inverse_kinematics (struct SixRSS * , const Ref<const Vector3f>& , const Ref<const Matrix3f>& , Ref<VectorXf>  );
void kinematic_inverse_jacobian (struct SixRSS * , const Ref<const Vector3f> & , const Ref<const Matrix3f> &, Ref<MatrixXf> & );
bool check_pose_and_compute_jacobian (struct SixRSS * , const Vector3f& , const Matrix3f& , Ref<MatrixXf> );
bool check_pose (struct SixRSS *, const Ref<const Vector3f> &, const Ref<const Matrix3f> &);
void generate_manipulator (struct SixRSS *);
void evaluate_workspace_discretization (struct SixRSS * , float *, float * );
void evaluate_workspace (int , float * , float * , float * );

#endif // MY_HEADER_H