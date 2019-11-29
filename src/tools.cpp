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
    std::cout << "ERROR - CalculateJacobian () - The state vector must have size 4." << "\n\n";
    return Hj;
  }
  
  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // TODO: YOUR CODE HERE 

  // pre-compute a set of terms to avoid repeated calculation
  double c1 = px*px+py*py;
  double c2 = sqrt(c1);
  double c3 = (c1*c2);

  // check division by zero
  if (fabs(c1) < 0.0001) {
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return Hj;
  }

  // compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
      -(py/c1), (px/c1), 0, 0,
      py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  std::cout << "tools.cpp: Calculated jacobian: \n" << x_state << "\n -> \n" << Hj << "\n\n";
  return Hj;
}

// This below is used in kalman_filter.cpp
VectorXd Tools::ConvertFromCartesianToPolar(const VectorXd& cart_coords) {
  
  // Extract cartesian coordinates from x
  double px = cart_coords(0);
  double py = cart_coords(1);
  double vx = cart_coords(2);
  double vy = cart_coords(3);

  if (px < 0.0000001) px = 0.0000001;
  if (py < 0.0000001) py = 0.0000001;

  double px2 = px * px;
  double py2 = py * py;
  
  // Magnitude
  double rho = sqrt(px2 + py2);
  if (rho < 0.0000001) rho = 0.0000001;
  
  //Angle
  double phi = atan2(py, px);
  
  // Rho dot
  double rhod = (px * vx + py * vy * 1.0) / rho;
  
  // Return results
  //VectorXd coords_polar(4);
  //coords_polar << rho, phi, rhod, 0;
  VectorXd coords_polar(3);
  coords_polar << rho, phi, rhod;

  std::cout << "tools.cpp: convert cart to polar: \n" << cart_coords << "\n -> \n" << coords_polar << "\n\n";
  return coords_polar;
}

// This below is used in Fu`sionEKF.cpp
VectorXd Tools::ConvertFromPolarToCartesian(const VectorXd& polar_coords) {
  
  // Extract polar coordinates
  // The Input file format is:
  // #R(radar) - meas_rho - meas_phi - meas_rho_dot - timestamp - gt_px - gt_py - gt_vx - gt_vy
  double rho = polar_coords(0);
  double phi = polar_coords(1);
  double rhod = polar_coords(2);
  
  // Clamp rho
  if (rho < 0.000001) rho = 0.000001;

  // Tranform from polar to cartesian
  double px = rho * cos(phi);
  double py = rho * sin(phi);
    
  // Output cartesian coordinates - no velocity
  VectorXd coords_cart(4);
  coords_cart << px, py, 0, 0;
  
  std::cout << "tools.cpp: convert polar to cart: \n" << polar_coords << " to " << coords_cart << "\n\n";
  return coords_cart;

}
