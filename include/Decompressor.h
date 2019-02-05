/*
  Decompressor code for the 2018 SCA paper "Collision-Aware and Online Compression of Rigid Body Simulations via Integrated Error Minimization"

  Timothy Jeruzalski - 2018
*/

#ifndef DECOMPRESSOR_H
#define DECOMPRESSOR_H

#include "CompressionEngine.h"

class Decompressor
{
public:
  Decompressor();

  // give a matrix of the compressed states for internal use
  bool setCompressedData(Eigen::MatrixXd* rotationData, Eigen::MatrixXd* positionData);

  CompressionEngine::Frame stateAtTime(double time);
private:
  Eigen::MatrixXd* rData;
  Eigen::MatrixXd* pData;

  // perform a binary search for the relevant time
  int32_t findIndexFromTime(const Eigen::MatrixXd* data, double time);
  // flip both quaternions into the same half sphere
  void flipQuats(Eigen::Quaterniond& q1, Eigen::Quaterniond& q2);

  Eigen::Quaterniond completeQuat(const Eigen::Vector4d& q);

  Eigen::Quaterniond decodeRotation(double time, int32_t index);
  Eigen::Vector3d    decodePosition(double time, int32_t index);
  
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif//DECOMPRESSOR_H
