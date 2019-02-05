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
  
  // because Eigen has "special" rotation behaviour:
  struct AngleAxis
  {
    Eigen::Vector3d ax;
    double an;
    
    AngleAxis()
    {
      ax << 1,0,0;
      an = 0;
    }
    AngleAxis(const Eigen::Quaternion<double>& q)
    {
      double a = 1 - q.w()*q.w();
      if(a != 0.0)
      {
	an = 2*atan2(sqrt(a), abs(q.w()));//q.w() >= 0 ? 2 * acos(q.w()) : 2.0*M_PI - 2*acos(q.w());
	ax << q.x() / sqrt(a), q.y() / sqrt(a), q.z() / sqrt(a);
	if(q.w() < 0)
	  ax = -ax;
      }
      else
      {
	an = 0;
	ax << 1, 0, 0;
      }      
    }
    double angle() const
    {
      return an;
    }
    double& angle()
    {
      return an;
    }
    Eigen::Vector3d axis() const
    {
      return ax;
    }
    Eigen::Vector3d& axis()
    {
      return ax;
    }
    Eigen::Quaternion<double> toQuat()
    {
      double ha = an / 2.0;
      return Eigen::Quaternion<double>( cos(ha), ax.x() * sin(ha), ax.y() * sin(ha), ax.z() * sin(ha) ).normalized();
    }
  };


private:
  // Very long method arguments -.-
  bool periodicitySearch(const std::vector<CompressionEngine::Frame>& intermediateFrames, const Eigen::Quaterniond& qFinal, double periodicity, const Eigen::AngleAxis<double>& AA, CompressionEngine::CompressedRotationFrame& fit);

  // since we only need the flattened E and W matrix for all of the operations, just need to construct it all flat
  bool constructFlatEMatrix(Eigen::VectorXd& E);
  bool constructFlatWMatrix(Eigen::VectorXd& W);
  // first construct the full W matrix, before flattening
  bool constructWMatrix(Eigen::Quaternion<double>& q);
  
  // the axis and angle for the full rotation between start and end frames
  Eigen::Quaternion<double> qInitial;
  Eigen::VectorXd eAxis;
  double eAngle;

  std::vector<Eigen::VectorXd> Emat;
  // the entries in the vector need to be scaled down by the number of elements in the vector when used
  std::vector<Eigen::VectorXd> Wmat;

  // the previously observed quaternion
  Eigen::Quaternion<double> PREV_QUAT;
  // reference to flip all of the quats around
  Eigen::Quaternion<double> REF_QUAT;  
};

#endif //ROTATIONINTERPOLATION_H
