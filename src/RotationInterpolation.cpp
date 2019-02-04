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
    Eigen::Quaterniond qMiddleLower = intermediateFrames[floor(intermediateFrames.size()/2.0)].rot;
    // df.getQuaternionFromState(&intermediateFrames[floor(intermediateFrames.size()/2.0) + 1], object, temp);
    Eigen::Quaterniond qMiddleUpper = intermediateFrames[floor(intermediateFrames.size()/2.0) + 1].rot;
    flipQuatToPrev(qMiddleLower);
    flipQuatToPrev(qMiddleUpper);
    drdf = 2*acos((qMiddleUpper * qMiddleLower.conjugate()).normalized().w());
  }
  double periodicity = drdf*intermediateFrames.size()/(2.0*M_PI);
  
  periodicitySearch(intermediateFrames, qFinal, periodicity, AA, fittedFrame);
  
  return fastLinearQuaternionErrorEstimate(intermediateFrames, fittedFrame);
}

Eigen::Quaternion<double> RotationInterpolation::completeQuatFromCompressedFrame(const Eigen::VectorXd& q)
{
  /*
    Takes a 3 component quat [x y z], and ups to a full quaternion [w x y z]
   */
  double wSq = (1 - q(0)*q(0) - q(1)*q(1) - q(2)*q(2));
  return Eigen::Quaternion<double>( wSq <= 0.001 ? 0 : sqrt(wSq) , q(0), q(1), q(2)).normalized();
}

double RotationInterpolation::fastLinearQuaternionErrorEstimate(const std::vector<CompressionEngine::Frame>& intermediateFrames, const CompressionEngine::CompressedRotationFrame& fittedFrame)
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
  PREV_QUAT.setIdentity();
  Eigen::Quaternion<double> qInitial, qFinal, qDiff;
  
  qFinal = completeQuatFromCompressedFrame(fittedFrame.rot.tail(3));
  
  qInitial = intermediateFrames.front().rot;


  flipQuat(qInitial);
  flipQuatToPrev(qFinal);

  qDiff = (qFinal * qInitial.conjugate()).normalized();

  Eigen::AngleAxis<double> AA(qDiff);
  
  for(uint32_t i = 0; i < checkIndices.size(); i++)
  {
    // get the normalized times:
    double t = checkIndices[i]/((double)intermediateFrames.size() - 1);

    AA.angle() = fittedFrame.rot(0) * t;

    // if the rotation is too small, the axis is wrong for numerics reasons (so do simple slerp in that case)
    Eigen::Quaterniond qInterp;
    if(abs(fittedFrame.rot(0)) > 0.1)
    {
      qInterp = Eigen::Quaternion<double>(AA) * qInitial;
    }
    else
    {
      qInterp = qInitial.slerp(t, qFinal);
	/*
	  (Eigen::VectorXd(4) << qInitial.x()*(1-t) + qFinal.x()*t,
	  qInitial.y()*(1-t) + qFinal.y()*t,
	  qInitial.z()*(1-t) + qFinal.z()*t,
	  qInitial.w()*(1-t) + qFinal.w()*t).finished().normalized();
	*/
    }
    
    Eigen::Quaterniond qState;

    qState = intermediateFrames[checkIndices[i]].rot;
    if(qState.coeffs().dot(qInterp.coeffs()) < 0)
    {
      qInterp.coeffs() *= -1.0;
    }
    error += ( qState.coeffs() - qInterp.coeffs() ).squaredNorm();
    
  }
  // normalize to as if all of the frames were considered
  return error * (intermediateFrames.size() / (double)checkIndices.size());       
}

bool RotationInterpolation::periodicitySearch(const std::vector<CompressionEngine::Frame>& intermediateFrames, const Eigen::Quaterniond& qFinal, double periodicity, const Eigen::AngleAxis<double>& AA, CompressionEngine::CompressedRotationFrame& fit)
{
  // do a search for the combinations of periodicities which gives the lowest cumulative error in the fit
  double bestError = -2;
  CompressionEngine::CompressedRotationFrame potentialFit(intermediateFrames.back().time);

  // to reconstruct the w from the norm identity, need w to be positive (ensured earlier)
  potentialFit.rot << 0, qFinal.x(), qFinal.y(), qFinal.z();
  
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
	potentialFit.rot(0) = -1.0*(2*M_PI - AA.angle() + std::floor(pOne)*2.0*M_PI);
      else
	potentialFit.rot(0) = (AA.angle() + std::floor(pOne)*2.0*M_PI);

      double error = fastLinearQuaternionErrorEstimate(intermediateFrames, potentialFit);

      if(error < bestError || bestError < -1)
      {
	fit = potentialFit;
	bestError = error;
      }
    }
  }
  return 1;
}

bool RotationInterpolation::flipQuat(Eigen::Quaternion<double>& q)
{
  if(REF_QUAT.dot(q) < 0)
  {
    q.coeffs() *= -1;
  }
  PREV_QUAT = q;
  return 1;
}

bool RotationInterpolation::flipQuatToPrev(Eigen::Quaternion<double>& q)
{
  if(PREV_QUAT.dot(q) < 0)
  {
    q.coeffs() *= -1;
  }
  PREV_QUAT = q;
  return 1;
}

