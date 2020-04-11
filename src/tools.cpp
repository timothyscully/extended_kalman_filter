#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if ((estimations.size() == 0) || (estimations[0].size() != ground_truth[0].size())) {
    cout << "CalculateRMSE() - Invalid size of estimations, or ground truth." << endl;
    return rmse;
  }

  for (int i = 0; i < estimations.size(); ++i) {
    VectorXd diff = estimations[i] - ground_truth[i];
    VectorXd sqdiff = diff.array() * diff.array();
    rmse += sqdiff;
  }

  rmse /= (float) (estimations.size());
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3, 4);

  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // TODO: YOUR CODE HERE
  // check division by zero
  if ((px == 0) && (py == 0)) {
    cout << "CalculateJacobian() - Divide by zero error" << endl;
    return Hj;
  }

  float px2py2 = px * px + py * py;
  float sqrtpx2py2 = sqrt(px2py2);
  float threetworootpx2py2 = powf(px2py2, 1.5f);

  // compute the Jacobian matrix
  Hj << (px / sqrtpx2py2), (py / sqrtpx2py2), 0, 0,
      -(py / px2py2), (px / px2py2), 0, 0,
      ((py * (vx * py - vy * px)) / threetworootpx2py2), ((px * (vy * px - vx * py)) / threetworootpx2py2),
      px / sqrtpx2py2, py / sqrtpx2py2;

  return Hj;
}
