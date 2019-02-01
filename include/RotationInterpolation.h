#ifndef ROTATIONINTERPOLATION_H
#define ROTATIONINTERPOLATION_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include "DataFormat.h"
#include <iostream>
#include <math.h>

class RotationInterpolation
{
public:
  RotationInterpolation();
  RotationInterpolation(DataFormat dataFormat, DataFormat compressedDataFormat)
  {
    PREV_QUAT = Eigen::Quaternion<double>(1,0,0,0);
    REF_QUAT = Eigen::Quaternion<double>(1,0,0,0);
    df = dataFormat;
    cDF = compressedDataFormat;
  }
  ~RotationInterpolation();

  double getError(double s);

  // computes the coefficients for the frame, and retuns the error DON"T USE
  double calculateFit(int32_t object, const std::vector<Eigen::VectorXd>& intermediateFrames, Eigen::VectorXd& fittedFrame);


  // computes the frame needed for the compression of the frames, and outputs in form [ α x y z ], where α is the angular velocity and [x y z] are the components from the quaternion
  double linearFit(int32_t object, const std::vector<Eigen::VectorXd>& intermediateFrames, Eigen::VectorXd& fittedFrame);
  // quadratic fitting has to use the other error metric (FULL_QUAT_ERROR())
  bool calculateQuadraticFit(int32_t object, const std::vector<Eigen::VectorXd>& intermediateFrames, Eigen::VectorXd& fittedFrame);

  bool interpolate(const Eigen::VectorXd& axisAngle, double normT, Eigen::VectorXd& quaternion);
  bool quadraticInterpolate(const Eigen::VectorXd& axisAngle, double normT, Eigen::VectorXd& quaternion);

  double fastQuadraticQuaternionErrorEstimate(const std::vector<Eigen::VectorXd>& intermediateFrames, const Eigen::Quaternion<double>& qFirst, const Eigen::MatrixXd& fittedFrame);
  double fastLinearQuaternionErrorEstimate(const std::vector<Eigen::VectorXd>& intermediateFrames, int32_t object, const Eigen::VectorXd& fittedFrame);  
  
  // double fastQuaternionErrorEstimate(const std::vector<Eigen::VectorXd>& intermediateFrames);
  bool axisAngleToQuat(const Eigen::VectorXd& axisAngle, double normT, Eigen::Quaternion<double>& q);
  Eigen::Quaternion<double> axisAngleToQuat(const Eigen::VectorXd& axisAngle, double normT);

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
  bool periodicitySearch(const std::vector<Eigen::VectorXd>& intermediateFrames, const Eigen::Quaternion<double>& qFirst, double periodicityOne, double periodicityTwo, const Eigen::MatrixXd& invT, const Eigen::AngleAxis<double>& AAone, const Eigen::AngleAxis<double>& AAtwo, Eigen::MatrixXd& fit);
  bool periodicitySearch(const std::vector<Eigen::VectorXd>& intermediateFrames, const Eigen::Quaternion<double>& qFinal, int32_t object, double periodicity, const Eigen::AngleAxis<double>& AA, Eigen::VectorXd& fit);
    
  bool quadraticAxisAngleToQuat(const Eigen::VectorXd& axisAngle, double normT, Eigen::Quaternion<double>& q);
  bool quadraticAxisAngleToQuat(const Eigen::VectorXd& wOne, const Eigen::VectorXd& wTwo, double normT, Eigen::Quaternion<double>& q);

  bool normSlerp(const Eigen::VectorXd& qOne, const Eigen::VectorXd& qTwo, Eigen::Quaternion<double>& quat);

  // determines the rotation between the prevQuat and the currently considered frame:
  //bool deltaRotation(const Eigen::Quaternion<double>& q);

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

  DataFormat df; // Is the uncompressed format in order to extract the quaternions from each state
  DataFormat cDF;
};

#endif //ROTATIONINTERPOLATION_H
