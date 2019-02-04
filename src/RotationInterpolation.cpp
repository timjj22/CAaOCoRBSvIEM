#include "RotationInterpolation.h"


RotationInterpolation::~RotationInterpolation()
{

}

double RotationInterpolation::getError(double s)
{
  if(Wmat.size() < 2)
    return 0;
  
  Eigen::VectorXd errorTerm, E, W;
  constructFlatEMatrix(E);
  constructFlatWMatrix(W);
  
  errorTerm = W - s*E;
  //printf(" ERROR: %f ", errorTerm.norm());
  return (errorTerm).norm();
}

// need to set up with a new frame from resetFrames() before calling
bool RotationInterpolation::constructWMatrix(Eigen::Quaternion<double>& q)
{
  // put the quaternion into the correct half sphere, and then put the axis angle into the Wmat
  flipQuat(q);

  // get the angle and axis for the difference between first and last frame
  Eigen::AngleAxis<double> qDiff(q * qInitial.conjugate());
  // to use Wmat need to scale by the size of the mat first
  Wmat.push_back(qDiff.angle() * qDiff.axis());

  eAxis  = qDiff.axis();
  eAngle = qDiff.angle();
  
  return 1;
}

bool RotationInterpolation::constructFlatEMatrix(Eigen::VectorXd& E)
{
  E = Eigen::VectorXd(Wmat.size() * 3);
  E.segment(0*Wmat.size(), Wmat.size()) = eAxis.x() * Eigen::VectorXd::Ones(Wmat.size());
  E.segment(1*Wmat.size(), Wmat.size()) = eAxis.y() * Eigen::VectorXd::Ones(Wmat.size());
  E.segment(2*Wmat.size(), Wmat.size()) = eAxis.z() * Eigen::VectorXd::Ones(Wmat.size());
  return 1;
}

bool RotationInterpolation::constructFlatWMatrix(Eigen::VectorXd& W)
{
  W = Eigen::VectorXd(Wmat.size() * 3);
  for(uint32_t i = 0; i < Wmat.size(); i++)
  {
    W(i + 0*Wmat.size()) = Wmat.size() * Wmat[i](0); 
    W(i + 1*Wmat.size()) = Wmat.size() * Wmat[i](1);
    W(i + 2*Wmat.size()) = Wmat.size() * Wmat[i](2);
  }
  return 1;
}

double RotationInterpolation::linearFit(const std::vector<CompressionEngine::Frame>& intermediateFrames, CompressionEngine::CompressedRotationFrame& fittedFrame)
{
  // get the angle and axis from the difference in quaternions from the first and last frame
  PREV_QUAT.setIdentity();

  Eigen::Quaterniond qFinal = intermediateFrames.back().rot;  
  flipQuat(qFinal);

  Eigen::Quaterniond qInitial = intermediateFrames.front().rot;
  flipQuat(qInitial);

  Eigen::Quaterniond qDiff = (qFinal * qInitial.conjugate()).normalized();
  Eigen::AngleAxis<double> AA(qDiff);

  double drdf = 0.0;
  if(intermediateFrames.size() > 3)
  {
    // get the periodicity
    // df.getQuaternionFromState(&intermediateFrames[floor(intermediateFrames.size()/2.0)], object, temp);    
    Eigen::Quaterniond qMiddleLower = intermediateFrames[floor(intermediateFrames.size()/2.0)];
    // df.getQuaternionFromState(&intermediateFrames[floor(intermediateFrames.size()/2.0) + 1], object, temp);
    Eigen::Quaterniond qMiddleUpper = intermediateFrames[floor(intermediateFrames.size()/2.0) + 1];
    flipQuatToPrev(qMiddleLower);
    flipQuatToPrev(qMiddleUpper);
    drdf = 2*acos((qMiddleUpper * qMiddleLower.conjugate()).normalized().w());
  }
  double periodicity = drdf*intermediateFrames.size()/(2.0*M_PI);
  
  periodicitySearch(intermediateFrames, qFinal, object, periodicity, AA, fittedFrame);
  
  return fastLinearQuaternionErrorEstimate(intermediateFrames, object, fittedFrame);
}

Eigen::Quaternion<double> RotationInterpolation::completeQuatFromCompressedFrame(const Eigen::VectorXd& q)
{
  /*
    Takes a 3 component quat [x y z], and ups to a full quaternion [w x y z]
   */
  double wSq = (1 - q(0)*q(0) - q(1)*q(1) - q(2)*q(2));
  return Eigen::Quaternion<double>( wSq <= 0.001 ? 0 : sqrt(wSq) , q(0), q(1), q(2)).normalized();
}

double RotationInterpolation::fastLinearQuaternionErrorEstimate(const std::vector<Eigen::VectorXd>& intermediateFrames, int32_t object, const Eigen::VectorXd& fittedFrame)
{
  double error = 0;
  
  std::vector<int32_t> checkIndices = {
    (int32_t)(0.07 * (intermediateFrames.size() - 1)),
    (int32_t)(0.22 * (intermediateFrames.size() - 1)),
    (int32_t)(0.55 * (intermediateFrames.size() - 1)),
    (int32_t)(0.63 * (intermediateFrames.size() - 1)),
    (int32_t)(0.72 * (intermediateFrames.size() - 1)),
    (int32_t)(0.89 * (intermediateFrames.size() - 1))
  };

  // get the axis for rotation
  //PREV_QUAT.setIdentity();
  Eigen::VectorXd temp;
  Eigen::Quaternion<double> qInitial, qFinal, qDiff;
  
  //df.getQuaternionFromState(&intermediateFrames.back(), object, temp);
  qFinal = completeQuatFromCompressedFrame(fittedFrame.tail(3));//Eigen::Quaternion<double>(temp(3), temp(0), temp(1), temp(2));
  //flipQuat(qFinal);

  df.getQuaternionFromState(&intermediateFrames[0], object, temp);
  qInitial = Eigen::Quaternion<double>(temp(3), temp(0), temp(1), temp(2));
  //flipQuat(qInitial);

  flipQuat(qInitial);
  flipQuatToPrev(qFinal);

  qDiff = (qFinal * qInitial.conjugate()).normalized();

  Eigen::AngleAxis<double> AA(qDiff);
  
  for(uint32_t i = 0; i < checkIndices.size(); i++)
  {
    // get the normalized times:
    double t = checkIndices[i]/((double)intermediateFrames.size() - 1);

    AA.angle() = fittedFrame(0) * t;

    // if the rotation is too small, the axis is wrong for numerics reasons (so do simple slerp in that case)
    Eigen::VectorXd qInterp;
    if(abs(fittedFrame(0)) > 0.1)
    {
      Eigen::Quaternion<double> interpQuaternion = Eigen::Quaternion<double>(AA) * qInitial;
      
      // member that XYZW
      qInterp = (Eigen::VectorXd(4) << interpQuaternion.x(), interpQuaternion.y(), interpQuaternion.z(), interpQuaternion.w()).finished();
    }
    else
    {
      qInterp = (Eigen::VectorXd(4) << qInitial.x()*(1-t) + qFinal.x()*t,
		 qInitial.y()*(1-t) + qFinal.y()*t,
		 qInitial.z()*(1-t) + qFinal.z()*t,
		 qInitial.w()*(1-t) + qFinal.w()*t).finished().normalized();
    }
    
    Eigen::VectorXd qState;
    for(int32_t j = 0; j < cDF.getNumberOfObjects(); j++)
    {
      df.getQuaternionFromState(&intermediateFrames[checkIndices[i]], j, qState);
      if(qState.dot(qInterp) < 0)
      {
	qInterp *= -1.0;
      }
      error += ( qState - qInterp ).squaredNorm();
    }
  }
  // normalize to as if all of the frames were considered
  return error * (intermediateFrames.size() / (double)checkIndices.size());       
}

bool RotationInterpolation::periodicitySearch(const std::vector<Eigen::VectorXd>& intermediateFrames, const Eigen::Quaternion<double>& qFinal, int32_t object, double periodicity, const Eigen::AngleAxis<double>& AA, Eigen::VectorXd& fit)
{
  // do a search for the combinations of periodicities which gives the lowest cumulative error in the fit
  double bestError = -2;
  Eigen::VectorXd potentialFit(4);

  // to reconstruct the w from the norm identity, need w to be positive (ensured earlier)
  potentialFit << 0, qFinal.x(), qFinal.y(), qFinal.z();
  
  for(int8_t j = 0; j <= 4; j++)
  {
    // do the search within a couple of periodicity cycles from the finite diff'ed values
    // start with offset of 0 though
    double offset = 0;
    if(intermediateFrames.size() > 5)
    {
      if(j == 1)
	offset = 1;
      else if(j == 2)
	offset = -1;
      else if(j == 3)
	offset = 2;
      else if(j == 4)
	offset = -2;
    }
    double pOne = periodicity + offset;

    for(uint8_t i = 0; i < 2; i++)
    {
      if(i == 1)
	potentialFit(0) = -1.0*(2*M_PI - AA.angle() + std::floor(pOne)*2.0*M_PI);
      else
	potentialFit(0) = (AA.angle() + std::floor(pOne)*2.0*M_PI);

      double error = fastLinearQuaternionErrorEstimate(intermediateFrames, object, potentialFit);
      //printf("error: %f\n", error);
      // penalize larger rotations for spins with large spacing between frames      
      //if(error < bestError || (error <= (bestError + eps) && abs(potentialFit(0)) < abs(fit(0)))  || bestError < -1)
      if(error < bestError || bestError < -1)
      {
	fit = potentialFit;
	bestError = error;
      }
    }
  }
  // if(bestError > 0.01)
  //   printf("  periodicity search best error: %f   ", bestError);
  return 1;
}

bool RotationInterpolation::normSlerp(const Eigen::VectorXd& qOne, const Eigen::VectorXd& qTwo, Eigen::Quaternion<double>& quat)
{
  Eigen::Quaternion<double> one(qOne(3), qOne(0), qOne(1), qOne(2));
  Eigen::Quaternion<double> two(qTwo(3), qTwo(0), qTwo(1), qTwo(2));

  flipQuat(one);
  flipQuat(two);

  quat = one.slerp(0.5, two);
  
  return 1;
}

// augmentedAxisAngle is the compressed axisangle and the initial quaternion stacked together: [angle, Ax, Ay, Az, Qx, Qy, Qz, Qw]
bool RotationInterpolation::interpolate(const Eigen::VectorXd& augmentedAxisAngle, double normT, Eigen::VectorXd& quaternion)
{
  // construct the quaternion from the axis angle:
  Eigen::Quaternion<double> qInterp;
  axisAngleToQuat(augmentedAxisAngle, normT, qInterp);
  
  // now rotate the quat WRT to the initial rotation
  Eigen::Quaternion<double> qFirst(augmentedAxisAngle(7), augmentedAxisAngle(4), augmentedAxisAngle(5), augmentedAxisAngle(6));

  qInterp = qInterp * qFirst;

  quaternion = (Eigen::VectorXd(4) << qInterp.x(), qInterp.y(), qInterp.z(), qInterp.w()).finished();

  return 1;
}

// axis angle is [alpha X Y Z]
Eigen::Quaternion<double> RotationInterpolation::axisAngleToQuat(const Eigen::VectorXd& axisAngle, double normT)
{
  return Eigen::Quaternion<double>(cos(normT*axisAngle(0)/2.0), axisAngle(1)*sin(normT*axisAngle(0)/2.0), axisAngle(2)*sin(normT*axisAngle(0)/2.0), axisAngle(3)*sin(normT*axisAngle(0)/2.0));
}
bool RotationInterpolation::axisAngleToQuat(const Eigen::VectorXd& axisAngle, double normT, Eigen::Quaternion<double>& q)
{
  q = Eigen::Quaternion<double>(cos(normT*axisAngle(0)/2.0), axisAngle(1)*sin(normT*axisAngle(0)/2.0), axisAngle(2)*sin(normT*axisAngle(0)/2.0), axisAngle(3)*sin(normT*axisAngle(0)/2.0));
  return 1;
}

bool RotationInterpolation::flipQuat(Eigen::Quaternion<double>& q)
{
  if(REF_QUAT.dot(q) < 0)
  {
    // there's no "-" overloaded for the quaternion operator 
    q.w() = -q.w();
    q.x() = -q.x();
    q.y() = -q.y();
    q.z() = -q.z();    
  }
  PREV_QUAT = q;
  return 1;
}

bool RotationInterpolation::flipQuatToPrev(Eigen::Quaternion<double>& q)
{
  if(PREV_QUAT.dot(q) < 0)
  {
    // there's no "-" overloaded for the quaternion operator 
    q.w() = -q.w();
    q.x() = -q.x();
    q.y() = -q.y();
    q.z() = -q.z();    
  }
  PREV_QUAT = q;
  return 1;
}

