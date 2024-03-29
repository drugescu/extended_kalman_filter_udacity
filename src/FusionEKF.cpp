#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09,   0, 0,
              0, 0.0009, 0,
              0, 0,   0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
    
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0,    0,    0,    // initial px known
             0, 1,    0,    0,    // initial py know
             0, 0, 1000,    0,    // initial vx unknown
             0, 0,    0, 1000;    // initial vy unknown

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << 0, 0, 0, 0, // dt4/4*noise_ax and all other elements are zero at first.
             0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0;

  H_laser_ << 1, 0, 0, 0,  // Select px and py for laser
              0, 1, 0, 0;

  // Initial state transition matrix F
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
  
  ekf_.H_ = H_laser_;
                
  std::cout << "Executed FusionEKF constructor, set up matrices." << std::endl;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  //std::cout << "ProcessMeasurement measurement_pack: \n" << measurement_pack.raw_measurements_ << std::endl << std::endl;
  //std::cout << "ProcessMeasurement sensor type: \n" << measurement_pack.sensor_type_ << std::endl << std::endl;
  //std::cout << "ProcessMeasurement timestamp: \n" << measurement_pack.timestamp_ << std::endl << std::endl;
  //fflush(stdout);
  
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    
    // Get inputs into a special vector
    // The Input file format is:
    //               0        1        2            3         4     5     6
    // #L(for laser) meas_px  meas_py  timestamp    gt_px     gt_py gt_vx gt_vy
    // #R(for radar) meas_rho meas_phi meas_rho_dot timestamp gt_px gt_py gt_vx gt_vy

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      
      VectorXd inputs(3);
      // Perform conversion to cartesian
    
      inputs << measurement_pack.raw_measurements_(0), 
                measurement_pack.raw_measurements_(1), 
                measurement_pack.raw_measurements_(2);
      ekf_.x_ << tools.ConvertFromPolarToCartesian(inputs);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0; // speeds 0,0
    }
    
    // timestamp
    previous_timestamp_ = measurement_pack.timestamp_ ;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
    // 0. Get timestamps first
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;   // in seconds
  previous_timestamp_  = measurement_pack.timestamp_;
  //cout << "FusionEKF.cpp: not first: dt = " << dt << " as in " << dt * 10000000.0 << " microseconds" << endl << endl;
  
  double noise_ax = 9.0;
  double noise_ay = 9.0;
  
  // 1. Modify the F matrix so that the time is integrated
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt,  0,
             0, 1,  0, dt,
             0, 0,  1,  0,
             0, 0,  0,  1;
    
  // 2. Set the process covariance matrix Q
  auto dt2 =  dt * dt;
  auto dt3 =  dt * dt2;
  auto dt4 = dt2 * dt2;
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt4/4.*noise_ax,     0,        dt3/2.*noise_ax,        0,
               0,            dt4/4.*noise_ay,         0,      dt3/2.*noise_ay,
               dt3/2.*noise_ax,     0,          dt2*noise_ax,        0,
               0,            dt3/2.*noise_ay,         0,        dt2*noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;

    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  //fflush(stdout);

}