#include "RotationInterpolation.h"


RotationInterpolation::~RotationInterpolation()
{

}

double RotationInterpolation::linearFit(const std::vector<CompressionEngine::Frame>& intermediateFrames, CompressionEngine::CompressedRotationFrame& fittedFrame)
{
  // get the angle and axis from the difference in quaternions from the first and last frame
  PREV_QUAT.setIdentity();

  Eigen::Quaterniond qFinal = intermediateFrames.back().rot;
  if(qFinal.w() < 0)
    qFinal.coeffs() *= -1.0;   
  // flipQuat(qFinal);

  Eigen::Quaterniond qInitial = intermediateFrames.front().rot;
  if(qInitial.w() < 0)
    qInitial.coeffs() *= -1.0;
  // flipQuat(qInitial);

  Eigen::Quaterniond qDiff = (qFinal * qInitial.conjugate()).normalized();
  Eigen::AngleAxis<double> AA(qDiff);

  double drdf = 0.0;
  if(intermediateFrames.size() >= 3)
  {
    // get the periodicity
    Eigen::Quaterniond qMiddleLower = intermediateFrames[floor(intermediateFrames.size()/2.0)].rot;
    Eigen::Quaterniond qMiddleUpper = intermediateFrames[floor(intermediateFrames.size()/2.0) + 1].rot;
    flipQuatToPrev(qMiddleLower);
    flipQuatToPrev(qMiddleUpper);
    drdf = 2*std::acos((qMiddleUpper * qMiddleLower.conjugate()).normalized().w());
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
  // PREV_QUAT.setIdentity();
  Eigen::Quaternion<double> qInitial, qFinal, qDiff;
  
  qFinal = completeQuatFromCompressedFrame(fittedFrame.rot.tail(3));
  
  qInitial = intermediateFrames.front().rot;

  // flip the qInitial so that the orientations are with +ve w:
  /*
  if(qInitial.w() < 0)
    qInitial.coeffs() *= -1.0;
  */
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
    if(std::abs(fittedFrame.rot(0)) > 0.1)
    {
      qInterp = Eigen::Quaternion<double>(AA) * qInitial;
    }
    else
    {
      // qInterp = qInitial.slerp(t, qFinal);
      qInterp.coeffs() = (1 - t) * qInitial.coeffs() + t * qFinal.coeffs();
      qInterp.normalize();
      
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
	fit.rot = potentialFit.rot;
	bestError = error;
      }
    }
  }
  return 1;
}

bool RotationInterpolation::flipQuat(Eigen::Quaternion<double>& q)
{
  if(REF_QUAT.coeffs().dot(q.coeffs()) < 0)
  {
    q.coeffs() *= -1;
  }
  PREV_QUAT = q;
  return 1;
}

bool RotationInterpolation::flipQuatToPrev(Eigen::Quaternion<double>& q)
{
  if(PREV_QUAT.coeffs().dot(q.coeffs()) < 0)
  {
    q.coeffs() *= -1;
  }
  PREV_QUAT = q;
  return 1;
}

