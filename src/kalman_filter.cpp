#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  // KF
  // x' = Fx + u
  x_ = F_ * x_;  
  
  // P' = F*P*Ft + Q
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;

  std::cout << "Predict :\n x_ = " << x_ << std::endl << "P_ = " << P_ << std::endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  // y = z - H * x'
  VectorXd y = z - H_ * x_;
  
  // S = H*P'*Ht + R
  MatrixXd Ht = H_.transpose();
  MatrixXd S  = H_ * P_ * Ht + R_;
  
  // K = P'*Ht*S^(-1)
  MatrixXd Sinv = S.inverse();
  MatrixXd K = P_ * Ht * Sinv;
  
  // x = x' + K*y
  x_ = x_ + (K * y);
  
  // P = (I - K*H) * P'
  MatrixXd I = MatrixXd::Identity(4, 4);
  P_ = (I - K * H_) * P_;
  std::cout << "Update KF:\n x_ = " << x_ << std::endl << "P_ = " << P_ << std::endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  Tools tools;
  
  //  KF             EKF
  // ------------------------
  // predict
  // x' = Fx + u    x'=f(x,u), u = 0
  // P' = FPFt+Q    use Fj instead of F

  // update
  // y=z-Hx'        y=h(x')
  MatrixXd Hj = tools.CalculateJacobian(x_);
  VectorXd y = z - tools.ConvertFromCartesianToPolar(x_);
  
  // Clamp values to -M_PI : M_PI
  // normalize the angle between -pi to pi
  //while( y(1) >  M_PI ) y(1) -= 2 * M_PI;
  //while( y(1) < -M_PI ) y(1) += 2 * M_PI;
  while (y(1) > M_PI || y(1) < -M_PI)
    if (y(1) > M_PI) y(1) -= M_PI;
    else y(1) += M_PI;
  
  // S=HP'Ht + R    use Hj instead of H
  MatrixXd Hjt = Hj.transpose();
  MatrixXd S  = Hj * P_ * Hjt + R_;
  
  // K=P'HtS-1
  MatrixXd Sinv = S.inverse();
  MatrixXd K = P_ * Hjt * Sinv;
  
  // x=x'+Ky
  x_ = x_ + (K * y);
  
  // P=(I-KH)P'
  MatrixXd I = MatrixXd::Identity(4, 4);
  P_ = (I - K * Hj) * P_;
  
  std::cout << "Update EKF:\n x_ = " << x_ << std::endl << "P_ = " << P_ << std::endl;
}
