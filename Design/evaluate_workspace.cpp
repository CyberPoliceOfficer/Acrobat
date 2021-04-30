#include "evaluate_workspace.hpp"


float get_home_z (struct SixRSS * manipulator){
    float A = pow(manipulator->r_m*cos(manipulator->theta_m(0)) - manipulator->r_b*cos(manipulator->theta_b(0)) + manipulator->h*cos(manipulator->beta_k(0)),2);
    float B = pow(manipulator->r_m*sin(manipulator->theta_m(0)) - manipulator->r_b*sin(manipulator->theta_b(0)) - manipulator->h*sin(manipulator->beta_k(0)),2);
    float C = A + B;
    float D = pow(manipulator->d,2);
    if (D - C > 0)
        return sqrt(D - C);
    else
        return -1;
} 

void compute_metrics (Ref<MatrixXf> Jacobian, float * LTCI, float * LRCI, float * LTSI, float * LRSI){
    JacobiSVD<MatrixXf> Jt(Jacobian.topRows<3>());
    *LTCI = Jt.singularValues()(2)/Jt.singularValues()(0);
    JacobiSVD<MatrixXf> Jr(Jacobian.bottomRows<3>());
    *LRCI = Jr.singularValues()(2)/Jr.singularValues()(0);
    *LTSI = Jacobian.topRows<3>().lpNorm<Infinity>();
    *LRSI = Jacobian.bottomRows<3>().lpNorm<Infinity>();
}

bool check_pose_and_compute_inverse_kinematic_jacobian (struct SixRSS * manipulator, const Vector3f& T, const Matrix3f& R, Ref<MatrixXf> Jacobian){
    Vector3f i_k;
    Vector<float, 6> alphas;

    float f_k, e_k, g_k, t_k;

    // Kinematic constraints
    for (int k = 0; k < 6; k++){
        
        i_k = T + R*manipulator->m_k.col(k) - manipulator->b_k.col(k);

        f_k = (cos(manipulator->beta_k(k))*i_k(0) + sin(manipulator->beta_k(k))*i_k(1))*2*manipulator->h;

        e_k = (sin(manipulator->beta_k(k))*sin(manipulator->phi_k(k))*i_k(0) - cos(manipulator->beta_k(k))*sin(manipulator->phi_k(k))*i_k(1) + cos(manipulator->phi_k(k))*i_k(2))*2*manipulator->h;
        
        g_k = i_k.squaredNorm() - (manipulator->d*manipulator->d - manipulator->h*manipulator->h);

        t_k = sqrt(e_k*e_k + f_k*f_k);
        if(abs(g_k) > t_k){
            return false;
        }

        alphas(k) = asin(g_k/t_k) - atan2(f_k, e_k);
    }

    // Joint constraints
    Matrix<float, 3, 6> h_k;
    Matrix<float, 3, 6> MR_k;
    Vector<float, 3> d_k;
    Vector<float, 3> j1_k;

    get_h_k(manipulator, alphas, h_k);
    MR_k = R*manipulator->m_k;

    for (int k = 0; k < 6; k++){
        d_k = MR_k.col(k) + T - h_k.col(k) - manipulator->b_k.col(k);
        j1_k = (h_k.col(k)).cross(d_k);
        float ang_j1_k = acos(j1_k.dot(manipulator->u_k.col(k))/(manipulator->u_k.col(k).norm()*j1_k.norm()));
        if (ang_j1_k > manipulator->joint_limit){
            return false;
        }
        float ang_j2_k = abs(acos(d_k.dot(MR_k.col(k))/(MR_k.col(k).norm()*d_k.norm())) - M_PI/2);
        if (ang_j2_k > manipulator->joint_limit){
            return false; 
        }
        float j_q = (h_k.col(k).cross(d_k)).dot(manipulator->u_k.col(k));
        Jacobian.row(k).head<3>() = d_k/j_q;
        Jacobian.row(k).tail<3>() = MR_k.col(k).cross(d_k)/j_q;
    }

    return true;
}

void get_h_k (struct SixRSS * manipulator, const Ref<const VectorXf>& alphas, Ref<MatrixXf> h_k){
    for (int k = 0; k < 6; k++){
        h_k(0,k) = manipulator->h*(sin(manipulator->beta_k(k))*sin(manipulator->phi_k(k))*sin(alphas(k)) + cos(manipulator->beta_k(k))*cos(alphas(k)));
        h_k(1,k) = manipulator->h*(-cos(manipulator->beta_k(k))*sin(manipulator->phi_k(k))*sin(alphas(k)) + sin(manipulator->beta_k(k))*cos(alphas(k)));
        h_k(2,k) = manipulator->h*(cos(manipulator->phi_k(k))*sin(alphas(k))); 
    }   
}

void inverse_kinematic_jacobian (struct SixRSS * manipulator, const Ref<const Vector3f>& T, const Ref<const Matrix3f>& R, Ref<MatrixXf> Jacobian){
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

    return;
}
 
bool check_pose (struct SixRSS * manipulator, const Vector3f& T, const Matrix3f& R){
    Vector3f i_k;
    Vector<float, 6> alphas;

    float f_k, e_k, g_k, t_k;

    // Kinematic constraints
    for (int k = 0; k < 6; k++){
        
        i_k = T + R*manipulator->m_k.col(k) - manipulator->b_k.col(k);

        f_k = (cos(manipulator->beta_k(k))*i_k(0) + sin(manipulator->beta_k(k))*i_k(1))*2*manipulator->h;

        e_k = (sin(manipulator->beta_k(k))*sin(manipulator->phi_k(k))*i_k(0) - cos(manipulator->beta_k(k))*sin(manipulator->phi_k(k))*i_k(1) + cos(manipulator->phi_k(k))*i_k(2))*2*manipulator->h;
        
        g_k = i_k.squaredNorm() - (manipulator->d*manipulator->d - manipulator->h*manipulator->h);

        t_k = sqrt(e_k*e_k + f_k*f_k);
        if(abs(g_k) > t_k){
            return false;
        }

        alphas(k) = asin(g_k/t_k) - atan2(f_k, e_k);
    }

    // Joint constraints
    Matrix<float, 3, 6> h_k;
    Matrix<float, 3, 6> MR_k;
    Vector<float, 3> d_k;
    Vector<float, 3> j1_k;

    get_h_k(manipulator, alphas, h_k);
    MR_k = R*manipulator->m_k;

    for (int k = 0; k < 6; k++){
        d_k = MR_k.col(k) + T - h_k.col(k) - manipulator->b_k.col(k);
        j1_k = (h_k.col(k)).cross(d_k);
        float ang_j1_k = acos(j1_k.dot(manipulator->u_k.col(k))/(manipulator->u_k.col(k).norm()*j1_k.norm()));
        if (ang_j1_k > manipulator->joint_limit){
            return false;
        }
        float ang_j2_k = abs(acos(d_k.dot(MR_k.col(k))/(MR_k.col(k).norm()*d_k.norm())) - M_PI/2);
        if (ang_j2_k > manipulator->joint_limit){
           return false; 
        }
    }

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
    
    float cell_volume = dx*dy*dz;
    int n_nodes = 0;
    float GTCI = 0, GRCI = 0, GTSI = 0, GRSI = 0;
    
    #pragma omp parallel for schedule(dynamic, 2) firstprivate(T) firstprivate(R) firstprivate(Jacobian) firstprivate(JacobianInv) reduction(+:GTCI) reduction(+:GRCI) reduction(+:GTSI) reduction(+:GRSI) reduction(+:n_nodes)
    for (int x = 0; x < n_x ; x++){
        T(0) = -1*xy_maj + x*dx;
        for (int y = 0; y < n_y ; y++){
            T(1) = -1*xy_maj + y*dy;
            for (int z = 0; z < n_z ; z++){                                
                T(2) = z_min + z*dz;
                if (check_pose_and_compute_inverse_kinematic_jacobian (evaluatee, T, R, JacobianInv)){                   
                    Jacobian = JacobianInv.inverse();
                    float LTCI, LRCI, LTSI, LRSI;
                    compute_metrics(Jacobian, &LTCI , &LRCI , &LTSI, &LRSI);
                    GTCI += LTCI;
                    GRCI += LRCI;
                    GTSI += LTSI;
                    GRSI += LRSI;
                    n_nodes++;
                }
            }
        }
    }                
    
    result[0] = n_nodes*cell_volume;
    result[1] = GTCI/n_nodes;
    result[2] = GRCI/n_nodes;
    result[3] = -1*GTSI/n_nodes;
    result[4] = -1*GRSI/n_nodes;
    return;
}
 
extern "C" {
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

    evaluatee.joint_limit = config[4];

    generate_manipulator(&evaluatee);

    switch (technique){ // technique 0 -> discretization, technique 1 -> tree
        case 0:
            evaluate_workspace_discretization (&evaluatee, config, result);
            break;
    } 

    return;
}
}
int main()
{
    float config[5] = {50, 50, 50, 3, M_PI/3};
    float params[8] = {0.5, 0.3, 0.2, 0.6667, 0.7, 0.3, 0.3491, M_PI/2}; //r_b, r_m, d_b, d_m, d, h, phi, beta
    float result[5] = {0, 0, 0, 0, 0}; //V, GTCI, GRCI, GTSI, GRSI
    evaluate_workspace(0, config, params, result);
    std::cout << "Volume = " << result[0] << std::endl;
    std::cout << "GTCI = " << result[1] << std::endl;
    std::cout << "GRCI = " << result[2] << std::endl;
    std::cout << "GTSI = " << result[3] << std::endl;
    std::cout << "GRSI = " << result[4] << std::endl;
}
