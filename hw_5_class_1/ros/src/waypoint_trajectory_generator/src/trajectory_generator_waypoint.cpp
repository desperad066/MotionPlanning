#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}
/*

    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function

    variable declaration: input       const int d_order,                    // the order of derivative
                                      const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
                                      const Eigen::MatrixXd &Vel,           // boundary velocity
                                      const Eigen::MatrixXd &Acc,           // boundary acceleration
                                      const Eigen::VectorXd &Time)          // time allocation in each segment
                          output      MatrixXd PolyCoeff(m, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients

*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::getQ(const Eigen::VectorXd &Time, const int p_order){
    int n_seg = Time.size();
    MatrixXd Q = MatrixXd::Zero((p_order+1)*n_seg, (p_order+1)*n_seg);
    
    for(int k = 0; k < n_seg; k++){
        MatrixXd Qk = MatrixXd::Zero(p_order+1, p_order+1);
        for(int i = 0; i<= p_order; i++){
            for(int j = 0; j <= p_order; j++){
                if(i < 4 || j < 4)
                    continue;
                else
                    Qk(i,j) = i*(i-1)*(i-2)*(i-3)*j*(j-1)*(j-2)*(j-3)*pow(Time(k), i+j-7)/(i+j-7);
            }
        }
        // Q.block<p_order+1, p_order+1>(k*(p_order+1), k*(p_order+1)) = Qk; whynot
        Q.block(k*(p_order+1), k*(p_order+1), p_order+1, p_order+1) = Qk;
    }
    return Q;

}
Eigen::MatrixXd TrajectoryGeneratorWaypoint::getM(const Eigen::VectorXd &Time, const int p_order){
    int n_seg = Time.size();
    MatrixXd M = MatrixXd::Zero(4*2*n_seg, (p_order+1)*n_seg);
    //4:p,v,a,j *2:start,end p_order+1:coeff
    for(int k = 0; k < n_seg; k++){
        MatrixXd Mk = MatrixXd::Zero(4*2, (p_order+1));
        Mk(0,0) = 1;
        Mk(1,1) = 1;
        Mk(2,2) = 2;
        Mk(3,3) = 6;

        double T = Time(k);
        // MatrixXd row1(1,p_order+1);
        // row1 << 1, T, pow(T,2), pow(T,3), pow(T,4), pow(T,5), pow(T,6), pow(T,7);
        // MatrixXd row2(1,p_order+1);
        // row2 << 0, 1, 2*T, 3*pow(T,2), 4*pow(T,3), 5*pow(T,4), 6*pow(T,5), 7*pow(T,6);
        // MatrixXd row3(1,p_order+1);
        // row3 << 0, 0, 2, 6*T, 12*pow(T,2), 20*pow(T,3), 30*pow(T,4), 42*pow(T,5);
        // MatrixXd row4(1,p_order+1);
        // row4 << 0, 0, 0, 6, 24*T, 60*pow(T,2), 120*pow(T,3), 210*pow(T,4);
        Mk.row(4) << 1, T, pow(T,2), pow(T,3), pow(T,4), pow(T,5), pow(T,6), pow(T,7);
        Mk.row(5) << 0, 1, 2*T, 3*pow(T,2), 4*pow(T,3), 5*pow(T,4), 6*pow(T,5), 7*pow(T,6);
        Mk.row(6) << 0, 0, 2, 6*T, 12*pow(T,2), 20*pow(T,3), 30*pow(T,4), 42*pow(T,5);
        Mk.row(7) << 0, 0, 0, 6, 24*T, 60*pow(T,2), 120*pow(T,3), 210*pow(T,4);

        M.block(k*(4*2), k*(p_order+1), 4*2, p_order+1) = Mk;
    }
    return M;


// function M = getM(n_seg, n_order, ts)
//     M = [];
//     for k = 1:n_seg
//         M_k = [];
//         %#####################################################
//         % STEP 1.1: calculate M_k of the k-th segment 
//         %
//         %
//         %
//         % M_k [4*2, 8] *2: start and end
//         M_k(1,:) = [1,0,0,0,0,0,0,0]; 
//         M_k(2,:) = [0,1,0,0,0,0,0,0];
//         M_k(3,:) = [0,0,2,0,0,0,0,0];
//         M_k(4,:) = [0,0,0,6,0,0,0,0];
        
//         T = ts(k);
//         M_k(5,:) = [1, T, T^2, T^3, T^4, T^5, T^6, T^7]; %p
//         M_k(6,:) = [0, 1, 2*T, 3*T^2, 4*T^3, 5*T^4, 6*T^5, 7*T^6]; %v
//         M_k(7,:) = [0, 0, 2, 6*T, 12*T^2, 20*T^3, 30*T^4, 42*T^5]; %a
//         M_k(8,:) = [0, 0, 0, 6, 24*T, 60*T^2, 120*T^3, 210*T^4]; %j
        
//         M = blkdiag(M, M_k);
//     end
// end

}


///
//choose matrix
///
Eigen::MatrixXd TrajectoryGeneratorWaypoint::getCt(const int n_seg, const int p_order){
    MatrixXd Ct = MatrixXd::Zero(n_seg*8, 8+4*(n_seg-1));

    MatrixXd eye = MatrixXd::Identity(4,4);

    Ct.block(0,0,4,4) = eye;
    for(int i = 1; i <= n_seg-1; i++){
        Ct(4+8*(i-1), 8+(i-1)) = 1;
        Ct.block(4+8*(i-1)+1, 8+n_seg-1+3*(i-1), 3, 3) = MatrixXd::Identity(3,3);

        Ct(4+8*(i-1)+4, 8+(i-1)) = 1;
        Ct.block(4+8*(i-1)+5, 8+n_seg-1+3*(i-1), 3, 3) = MatrixXd::Identity(3,3);
    }
    Ct.block(4+8*(n_seg-1), 4, 4, 4) = eye;
    return Ct;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial =====?????===== the difference of the order of polynomial and the order of derivative
    int p_num1d   = p_order + 1;                  // the number of variables in each segment

    int n_seg = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(n_seg, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * n_seg), Py(p_num1d * n_seg), Pz(p_num1d * n_seg);
    bool debug = false;
    ROS_ERROR("Time");
    if(debug) std::cout << Time << std::endl;

    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */
    ROS_ERROR("Q");
    Eigen::MatrixXd Q = getQ(Time, p_order);
    if(debug) std::cout << Q << std::endl;

    ROS_ERROR("M");
    Eigen::MatrixXd M = getM(Time, p_order);
    if(debug) std::cout << M << std::endl;

    ROS_ERROR("Ct");
    Eigen::MatrixXd Ct = getCt(n_seg, p_order);
    if(debug) std::cout << Ct << std::endl;

    Eigen::MatrixXd C = Ct.transpose();
    
    ROS_ERROR("R");
    Eigen::MatrixXd R = C*M.inverse().transpose()*Q*M.inverse()*Ct;
    if(debug) std::cout << R << std::endl;

    ROS_ERROR("R_pp");
    Eigen::MatrixXd R_pp = R.block(n_seg+7, n_seg+7, 3*(n_seg-1), 3*(n_seg-1));
    if(debug) std::cout << R_pp << std::endl;
    
    ROS_ERROR("R_fp");
    Eigen::MatrixXd R_fp = R.block(0, n_seg+7, n_seg+7, 3*(n_seg-1));
    if(debug) std::cout << R_fp << std::endl;

    Eigen::VectorXd dF[3];//X, Y ,Z
    Eigen::VectorXd dP[3];//X, Y ,Z
    Eigen::VectorXd dFplusdP[3];
    Eigen::VectorXd coeff[3];//X, Y ,Z
    
    ROS_ERROR("dF");
    for(int i = 0; i < 3; i++){
        dF[i] = VectorXd::Zero(n_seg+7);
        dF[i](0) = Path(0,i);
        dF[i](4) = Path(Path.rows()-1, i);
        
        for(int j =0; j< n_seg-1; j++)
            dF[i](8+j) = Path(1+j, i);
        dP[i] = -R_pp.inverse()* R_fp.transpose() * dF[i];
        dFplusdP[i] = VectorXd::Zero(4*(n_seg+1));
        dFplusdP[i].head(n_seg+7) = dF[i];
        dFplusdP[i].tail(3*(n_seg-1)) = dP[i];
        ROS_ERROR("coedd,%d, %d, %d", i, M.rows(), M.cols());
        if(debug) std::cout << dFplusdP[i] << std::endl;
        ROS_ERROR("coeff");
        coeff[i] = M.inverse()*Ct*dFplusdP[i];
        if(debug) std::cout << coeff[i] << std::endl;

    }

    ROS_ERROR("PolyCoeff");
    for(int i = 0; i < n_seg; i++){
        RowVectorXd PolyCoeff_seg(3 * p_num1d);
        for(int j = 0; j < p_num1d; j ++){
            PolyCoeff_seg(j) = coeff[0](i*p_num1d+j);
            PolyCoeff_seg(j+p_num1d) = coeff[1](i*p_num1d+j);
            PolyCoeff_seg(j+2*p_num1d) = coeff[2](i*p_num1d+j);
        }
        PolyCoeff.row(i) = PolyCoeff_seg;
    }
    if(debug) std::cout << PolyCoeff << std::endl;
    
    
    // =====Matlab=====
    // C = Ct';
    // R = C * inv(M)' * Q * inv(M) * Ct;
    // R_cell = mat2cell(R, [n_seg+7 3*(n_seg-1)], [n_seg+7 3*(n_seg-1)]);
    // R_pp = R_cell{2, 2};
    // R_fp = R_cell{1, 2};
    // dF = [start_cond';end_cond';waypoints(2:end-1)];
    // dP = -inv(R_pp) * R_fp' * dF;
    // poly_coef = inv(M) * Ct * [dF;dP];





    /*   Produce the dereivatives in X, Y and Z axis directly.  */
    












    /*   Produce the Minimum Snap cost function, the Hessian Matrix   */
    










    return PolyCoeff;
}
