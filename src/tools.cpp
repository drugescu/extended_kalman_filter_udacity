#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  
  if ( x_state.size() != 4 ) {
    std::cout << "ERROR - CalculateJacobian () - The state vector must have size 4." << std::endl;
    return Hj;
  }
  
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // TODO: YOUR CODE HERE 

  // pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  // check division by zero
  if (fabs(c1) < 0.0001) {
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return Hj;
  }

  // compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
      -(py/c1), (px/c1), 0, 0,
      py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  std::cout << "calculate jacobian: " << Hj << std::endl;
  return Hj;
}

VectorXd Tools::ConvertFromCartesianToPolar(const VectorXd& cart_coords) {
  
  // Extract cartesian coordinates from x
  auto px = cart_coords(0);
  auto py = cart_coords(1);
  auto vx = cart_coords(2);
  auto vy = cart_coords(3);
  
  auto px2 = px * px;
  auto py2 = py * py;
  
  // Magnitude
  auto rho = sqrt(px2 + py2);
  
  //Angle
  auto phi = atan2(py, px);
  
  // Rho dot
  auto rhod = (px * vx + py * vy * 1.0) / rho;
  
  // Return results
  auto coords_polar = VectorXd(3);
  coords_polar << rho, phi, rhod;

  std::cout << "convert cart to polar: " << cart_coords << "," << coords_polar << std::endl;
  return coords_polar;
}

VectorXd Tools::ConvertFromPolarToCartesian(const VectorXd& polar_coords) {
  
  // Extract polar coordinates
  // The Input file format is:
  // #R(radar) - meas_rho - meas_phi - meas_rho_dot - timestamp - gt_px - gt_py - gt_vx - gt_vy
  auto rho = polar_coords(0);
  auto phi = polar_coords(1);
  auto rhod = polar_coords(2);
  
  // Clamp rho
  if (rho < 0.00001) rho = 0.00001;

  // Tranform from polar to cartesian
  auto px = rho * cos(phi);
  auto py = rho * sin(phi);
    
  // Output cartesian coordinates - no velocity
  auto coords_cart = VectorXd(4);
  coords_cart << px, py, 0, 0;
  
  return coords_cart;

}