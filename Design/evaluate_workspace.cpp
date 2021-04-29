#include "evaluate_workspace.hpp"

void get_h_k (struct SixRSS * manipulator, const Ref<const VectorXf>& alphas, Ref<MatrixXf> h_k){
    for (int k = 0; k < 6; k++){
        h_k(0,k) = manipulator->h*(sin(manipulator->beta_k(k))*sin(manipulator->phi_k(k))*sin(alphas(k)) + cos(manipulator->beta_k(k))*cos(alphas(k)));
        h_k(1,k) = manipulator->h*(-cos(manipulator->beta_k(k))*sin(manipulator->phi_k(k))*sin(alphas(k)) + sin(manipulator->beta_k(k))*cos(alphas(k)));
        h_k(2,k) = manipulator->h*(cos(manipulator->phi_k(k))*sin(alphas(k))); 
    }   
}

void inverse_kinematic_jacobian (struct SixRSS * manipulator, const Ref<const Vector3f>& T, const Ref<const Matrix3f>& R, Ref<MatrixXf>  Jacobian){
    /*
    MATLAB:
    M_k = R*m_k + T;
    Jx = [(M_k-H_k)' cross(M_k-T, M_k-H_k)'];
    Jq = eye(6).*dot(cross(H_k-b_k,M_k - H_k),u_k);
    J = inv(Jx)*Jq;
    */
    Matrix<float, 3, 6> MR_k;
    Matrix<float, 3, 6> h_k;
    Vector<float, 6> alphas;
    Vector<float, 3> d_k;
    
    inverse_kinematics(manipulator, T, R, alphas);
    get_h_k(manipulator, alphas, h_k);
    MR_k = R*manipulator->m_k;
    
    for (int k = 0; k < 6; k++){
        d_k = MR_k.col(k) + T - h_k.col(k) - manipulator->b_k.col(k);
        float j_q = (h_k.col(k).cross(d_k)).dot(manipulator->u_k.col(k));
        Jacobian.row(k).head<3>() = d_k/j_q;
        Jacobian.row(k).tail<3>() = MR_k.col(k).cross(d_k)/j_q;
    }
}

void inverse_kinematics (struct SixRSS * manipulator, const Ref<const Vector3f>& T, const Ref<const Matrix3f>& R, Ref<VectorXf> alphas){

    Vector3f i_k;

    float f_k, e_k, g_k, t_k;
    
    // Kinematic constraints
    for (int k = 0; k < 6; k++){
        
        i_k = T + R*manipulator->m_k.col(k) - manipulator->b_k.col(k);

        f_k = (cos(manipulator->beta_k(k))*i_k(0) + sin(manipulator->beta_k(k))*i_k(1))*2*manipulator->h;

        e_k = (sin(manipulator->beta_k(k))*sin(manipulator->phi_k(k))*i_k(0) - cos(manipulator->beta_k(k))*sin(manipulator->phi_k(k))*i_k(1) + cos(manipulator->phi_k(k))*i_k(2))*2*manipulator->h;

        g_k = i_k.squaredNorm() - (manipulator->d*manipulator->d - manipulator->h*manipulator->h);

        t_k = sqrt(e_k*e_k + f_k*f_k);
        if(abs(g_k) <= t_k){
            alphas(k) = asin(g_k/t_k) - atan2(f_k, e_k);
        }
    }

}
 
bool check_pose (struct SixRSS * manipulator, const Vector3f& T, const Matrix3f& R){
    Vector3f i_k;

    float f_k, e_k, g_k;

    // Kinematic constraints
    for (int k = 0; k < 6; k++){
        
        i_k = T + R*manipulator->m_k.col(k) - manipulator->b_k.col(k);

        f_k = (cos(manipulator->beta_k(k))*i_k(0) + sin(manipulator->beta_k(k))*i_k(1))*2*manipulator->h;

        e_k = (sin(manipulator->beta_k(k))*sin(manipulator->phi_k(k))*i_k(0) - cos(manipulator->beta_k(k))*sin(manipulator->phi_k(k))*i_k(1) + cos(manipulator->phi_k(k))*i_k(2))*2*manipulator->h;
        
        g_k = i_k.squaredNorm() - (manipulator->d*manipulator->d - manipulator->h*manipulator->h);

        if(abs(g_k) > sqrt(e_k*e_k + f_k*f_k)){
            return false;
        }
    }

    // Joint constraints

    return true;
}

void generate_manipulator (struct SixRSS * manipulator){
    int k;
    int n;
    for (int i = 0; i < 6 ; i++){
        k = i + 1;
        n = floor(i/2);
        manipulator->theta_b(i) = n*(2.0/3)*M_PI + pow(-1,k)*asin(manipulator->d_b);
        manipulator->theta_m(i) = n*(2.0/3)*M_PI + pow(-1,k)*asin(manipulator->d_m);
        
        manipulator->b_k(0,i) = manipulator->r_b * cos(manipulator->theta_b(i));
        manipulator->b_k(1,i) = manipulator->r_b * sin(manipulator->theta_b(i));
        manipulator->b_k(2,i) = 0;

        manipulator->m_k(0,i) = manipulator->r_m * cos(manipulator->theta_m(i));
        manipulator->m_k(1,i) = manipulator->r_m * sin(manipulator->theta_m(i));
        manipulator->m_k(2,i) = 0;
            
        manipulator->beta_k(i) = n*(2.0/3)*M_PI + pow(-1,k)*manipulator->beta;
        manipulator->phi_k(i) = pow(-1,k+1)*manipulator->phi;

        manipulator->u_k(0,i) = -cos(manipulator->beta_k(i) + M_PI/2)*cos(manipulator->phi_k(i));
        manipulator->u_k(1,i) = -sin(manipulator->beta_k(i) + M_PI/2)*cos(manipulator->phi_k(i));
        manipulator->u_k(2,i) = -sin(manipulator->phi_k(i));       
    }

}

void evaluate_workspace_discretization (struct SixRSS * evaluatee, float * config, float * result){

    float xy_maj = evaluatee->d + evaluatee->h - cos(M_PI/3)*evaluatee->r_b;
    float z_maj = evaluatee->d + evaluatee->h;
    float z_min = 0;
    
    int n_x = static_cast<int>(config[0]);
    int n_y = static_cast<int>(config[1]);
    int n_z = static_cast<int>(config[2]);

    float dx = 2*xy_maj/(n_x-1);
    float dy = 2*xy_maj/(n_y-1);
    float dz = z_maj/(n_z-1);
    
    // Pose
    Vector3f T;   
    Matrix3f R; //R = AngleAxisf(pose(3), Vector3f::UnitZ()) * AngleAxisf(pose(4), Vector3f::UnitY()) * AngleAxisf(pose(5), Vector3f::UnitZ());
    R << 1, 0, 0,
         0, 1, 0,
         0, 0, 1;
    
    Matrix<float, 6, 6> Jacobian;
    Matrix<float, 6, 6> JacobianInv;
    
    T << 0, 0, 0.8;
    inverse_kinematic_jacobian (evaluatee, T, R, Jacobian);
    std::cout << Jacobian << std::endl;

    float cell_volume = dx*dy*dz;
    int n_nodes = 0;

    for (int x = 0; x < n_x ; x++){
        T(0) = -1*xy_maj + x*dx;
        for (int y = 0; y < n_y ; y++){
            T(1) = -1*xy_maj + y*dy;
            for (int z = 0; z < n_z ; z++){                                
                T(2) = z_min + z*dz;
                if (check_pose(evaluatee, T, R)){
                    inverse_kinematic_jacobian (evaluatee, T, R, Jacobian);
                    Jacobian = Jacobian.inverse();
                    n_nodes++;
                }
            }
        }
    }                
    
    std::cout << n_nodes << std::endl;
    result[0] = n_nodes*cell_volume;
    return;
}
 
void evaluate_workspace (int technique, float * config, float * params, float * result){

    struct SixRSS evaluatee;

    evaluatee.r_b = params[0]; //Radious of the base
    evaluatee.r_m = params[1]; //Radious of the platform
        
    evaluatee.d_b = params[2]; //Angle between base's anchors
    evaluatee.d_m = params[3]; //Angle between platform's anchors
        
    evaluatee.d = params[4]; //Platform's arm lenght
    evaluatee.h = params[5]; //Servo's arm lenght
        
    evaluatee.phi = params[6]; //Angle between servo's arm and platform's base
    evaluatee.beta = params[7];

    generate_manipulator(&evaluatee);

    switch (technique){ // technique 0 -> discretization, technique 1 -> tree
        case 0:
            evaluate_workspace_discretization (&evaluatee, config, result);
            break;
    } 

    return;
}


int main()
{
  float config[4] = {50, 50, 50, 3};
  float params[8] = {0.5, 0.3, 0.2, 0.6667, 0.7, 0.3, 0.3491, M_PI/2}; //r_b, r_m, d_b, d_m, d, h, phi, beta
  float result[5] = {0, 0, 0, 0, 0};
  evaluate_workspace(0, config, params, result);
  std::cout << result[0] << std::endl;
}
