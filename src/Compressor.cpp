/*
  Compressor code for the 2018 SCA paper "Collision-Aware and Online Compression of Rigid Body Simulations via Integrated Error Minimization"

  Timothy Jeruzalski - 2018
 */

#include "Compressor.h"

Compressor::Compressor(double positionThreshold/* = 5e-2*/,
		       double positionPeakThreshold/* = 5e-4*/,
		       double rotationThreshold/* = 1e-2*/,
		       double rotationPeakThreshold/* = 5e-3*/)
{
  pThreshold  = positionThreshold;
  pPThreshold = positionPeakThreshold;
  rThreshold  = rotationThreshold;
  rPThreshold = rotationPeakThreshold;

  qSolve = Quadratic();
}

Compressor::~Compressor()
{

}

bool Compressor::compressFrame(const Frame& newFrame, bool force /* = 0*/)
{
  intermediatePositionFrames.push_back(newFrame);
  intermediateRotationFrames.push_back(newFrame);

  return compressFramePosition(force) || compressFrameRotation(force);
}

bool Compressor::compressFramePosition(bool force)
{
  bool compressed = 0;
  // the position step:
  qSolve.fastPrecompute(intermediatePositionFrames);
  Eigen::MatrixXd weights;
  qSolve.fastSolve((Eigen::MatrixXd(3, 2) << intermediatePositionFrames.front().pos, intermediateFrames.back().pos).finished().transpose(), weights);
  positionError = qSolve.fastErrorEstimate(weights);
  
  // now set the current compressed frame
  cPosFrame = CompressionEngine::CompressedPositionFrame(intermediatePositionFrames.back().time, (Eigen::VectorXd(6) << weights.row(1).transpose(), intermediatePositionFrames.back().pos).finished());
  
  // keyframe placement if error bounds exceeded, or if forced for the current keyframe
  if(contactFlag || force || compressedPositionFrames.size() == 0)
  {
    // add the current frame
    compressedPositionFrames.push_back(cPosFrame);
    compressed = 1;

    // keep the current frame in the window
    intermediatePositionFrames.erase(intermediatePositionFrames.begin(), intermediatePositionFrames.end() - 1);

    // remove the contact flag:
    contactFlag = 0;
  }
  else if(peakDetected(pPThreshold, positionError, positionErrorPrev, positionErrorDiffPrev))
  {
    // upon a peak, add the f-2 keyframe:
    compressedPositionFrames.push_back(cppPosFrame);
    compressed = 1;
   
    // clear the old intermediate frames out and reset the error averages
    intermediatePositionFrames.erase(intermediatePositionFrames.begin(), intermediatePositionFrames.end() - 3);
    qSolve.resetFit(&intermediatePositionFrames);

    // now set the flag to drop the current frame next run through
    contactFlag = 1;

  }
  else if(std::abs(positionError) > pThreshold)
  {
    // upon error bounds tripping, add the f-1 keyframe:
    compressedPositionFrames.push_back(cpPosFrame);
    compressed = 1;

    //
  }
  
  // update the previous frames
  cppPosFrame = cpPosFrame;  
  cpPosFrame = cPosFrame;
  
  return compressed;
}

bool Compressor::compressFrameRotation(bool force)
{  
  // the rotation step:
  Eigen::VectorXd frame;
  rotationError = 0.0;
  compressed = 0;
  
  if(intermediateFrames.size() >= 2)
    rotationError = RI.linearFit(intermediateRotationFrames, frame);
  else
    frame = (Eigen::VectorXd(4) << 0, intermediateRotationFrames.back().rot.x(), intermediateRotationFrames.back().rot.y(), intermediateRotationFrames.back().rot.z()).finished();

  // now copy into the current frame
  cRotFrame = CompressionEngine::CompressedRotationFrame(intermediateRotationFrames.back().time, frame);

  // keyframe placement:
  if(force || compressedRotationFrames.size() == 0)
  {
    compressedRotationFrames.push_back(cRotFrame);
    compressed = 1;
  }
  else if(peakDetected(rPThreshold, rotationError, rotationErrorPrev, rotationErrorDiffPrev))
  {
    compressedRotationFrames.push_back(cppRotFrame);
    compressed = 1;
  }
  else if(std::abs(rotationError) > rThreshold)
  {
    compressedRotationFrames.push_back(cpRotFrame);
    compressed = 1;
  }

  cppRotFrame = cpRotFrame;
  cpRotFrame = cRotFrame;

  return compressed;
}

bool Compressor::peakDetected(const double& threshold, const double& error, double& prevError, double& prevErrorDiff)
{
  /*
    Performs a finite difference second derivative approximation on the error signal
    
    d²e/dt² ≈ ((eₜ - eₜ₋₁) - (eₜ₋₁ - eₜ₋₂)) / dt²
  */
  bool peak = 0;
  double currentDiff = error - prevError; // ≈ de/dt

  if(std::abs(currentDiff - prevErrorDiff) > threshold && intermediateFrames.size() > 2)
  {
    // then we have a peak, reset previous
    peak = 1;
    prevError = 0;
    prevErrorDiff = 0;
  }
  else
  {
    peak = 0;
    prevError = error;
    prevErrorDiff = currentDiff;
  }
  return peak;
}

Eigen::MatrixXd Compressor::getCompressedMatrixPositionRepresentation()
{
  Eigen::MatrixXd cMat(7, compressedPositionFrames.size());
  for(int32_t i = 0; i < CompressedPositionFrames.size(); i++)
  {
    cMat(0, i) = compressedPositionFrames[i].time;
    cMat.col(i).tail(6) = compressedPositionFrames[i].pos;
  }

  return cMat;
}

Eigen::MatrixXd Compressor::getCompressedMatrixRotationRepresentation()
{
  Eigen::MatrixXd cMat(5, compressedRotationFrames.size());
  for(int32_t i = 0; i < CompressedRotationFrames.size(); i++)
  {
    cMat(0, i) = compressedRotationFrames[i].time;
    cMat.col(i).tail(4) = compressedRotationFrames[i].rot;
  }

  return cMat;
}

std::vector<CompressionEngine::CompressedPositionFrame> Compressor::getCompressedVectorPositionRepresentation()
{
  return compressedPositionFrames;
}

std::vector<CompressionEngine::CompressedRotationFrame> Compressor::getCompressedVectorRotationRepresentation()
{
  return compressedRotationFrames;
}
