/*
  Rotation handling for the compression system

  Timothy Jeruzalski - 2018
*/

#ifndef ROTATIONINTERPOLATION_H
#define ROTATIONINTERPOLATION_H

#include "CompressionEngine.h"

class RotationInterpolation
{
public:
  RotationInterpolation()
  {
    PREV_QUAT = Eigen::Quaternion<double>(1,0,0,0);
    REF_QUAT = Eigen::Quaternion<double>(1,0,0,0);
  }
  ~RotationInterpolation();

  double getError(double s);
  
  // computes the frame needed for the compression of the frames, and outputs in form [ α x y z ], where α is the angular velocity and [x y z] are the components from the quaternion
  double linearFit(const std::vector<CompressionEngine::Frame>& intermediateFrames, CompressionEngine::CompressedRotationFrame& fittedFrame);
    
  double fastLinearQuaternionErrorEstimate(const std::vector<CompressionEngine::Frame>& intermediateFrames, const CompressionEngine::CompressedRotationFrame& fittedFrame);  

  Eigen::Quaternion<double> completeQuatFromCompressedFrame(const Eigen::VectorXd& q);

  // Compares the quaternion to [1 0 0 0] quaternion to determine if the quat needs to be inverted
  bool flipQuat(Eigen::Quaternion<double>& q);
  // compare to previously observed quaternion
  bool flipQuatToPrev(Eigen::Quaternion<double>& q);

private:
  // Very long method arguments -.-
  bool periodicitySearch(const std::vector<CompressionEngine::Frame>& intermediateFrames, const Eigen::Quaterniond& qFinal, double periodicity, const Eigen::AngleAxis<double>& AA, CompressionEngine::CompressedRotationFrame& fit);
  // the previously observed quaternion
  Eigen::Quaternion<double> PREV_QUAT;
  // reference to flip all of the quats around
  Eigen::Quaternion<double> REF_QUAT;  
};

#endif //ROTATIONINTERPOLATION_H
