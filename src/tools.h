#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"
#include <iostream>

class Tools {
 public:
  /**
   * Constructor.
   */
  Tools();

  /**
   * Destructor.
   */
  virtual ~Tools();

  /**
   * A helper method to calculate RMSE.
   */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, 
                                const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   * A helper method to calculate Jacobians.
   */
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

  /**
   * My own helper method to convert from Cartesian to Polar coordinates to calculate y for radar
   */
  Eigen::VectorXd ConvertFromCartesianToPolar(const Eigen::VectorXd& cart_coords);
  
  /**
   * My own helper method to convert from Polar to Cartesian coordinates to calculate actual z from radar data
   */
  Eigen::VectorXd ConvertFromPolarToCartesian(const Eigen::VectorXd& polar_coords);

  const double EPS = 0.0001;
};

#endif  // TOOLS_H_