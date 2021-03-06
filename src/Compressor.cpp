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
  RI = RotationInterpolation();
}

Compressor::~Compressor()
{

}

bool Compressor::compressFrame(const CompressionEngine::Frame& newFrame, bool contact, bool force /* = 0*/)
{
  intermediatePositionFrames.emplace_back(newFrame);
  intermediateRotationFrames.emplace_back(newFrame);

  // if we know there is a contact starting now, then drop this frame and the this+2 frame keyframes
  bool forceFrame = force;
  if(contact)
  {
    // force this frame, and the one 2 in the future
    futureKeyframe = 2;
    forceFrame = 1;
  }
  else
  {
    // decrement the future keyframe count, but don't worry about wrap-around
    futureKeyframe = std::max(-1, futureKeyframe - 1);
  }
  if(futureKeyframe == 0)
  {
    forceFrame |= 1;
  }

  return compressFramePosition(forceFrame) | compressFrameRotation(forceFrame); // the || short circuits 
}

bool Compressor::compressFramePosition(bool force)
{
  bool compressed = 0;
  // the position step:
  qSolve.fastPrecompute(intermediatePositionFrames);
  Eigen::MatrixXd weights;
  qSolve.fastSolve((Eigen::MatrixXd(3, 2) << intermediatePositionFrames.front().pos, intermediatePositionFrames.back().pos).finished().transpose(), weights);
  positionError = qSolve.fastErrorEstimate(weights);
  
  // now set the current compressed frame
  cPosFrame = CompressionEngine::CompressedPositionFrame(intermediatePositionFrames.back().time, (Eigen::VectorXd(6) << weights.row(1).transpose(), intermediatePositionFrames.back().pos).finished());
  
  // keyframe placement if error bounds exceeded, or if forced for the current keyframe
  if(force || compressedPositionFrames.size() == 0)
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
    printf("pos peak\n");
    // upon a peak, add the f-2 keyframe:
    compressedPositionFrames.push_back(cppPosFrame);
    compressed = 1;
   
    // clear the old intermediate frames out and reset the error averages
    intermediatePositionFrames.erase(intermediatePositionFrames.begin(), intermediatePositionFrames.end() - 3);
    Eigen::MatrixXd iFit = qSolve.resetFit(&intermediatePositionFrames);

    // add in the f-2 -> f keyframe
    /*
    qSolve.fastPrecompute(intermediatePositionFrames);
    qSolve.fastSolve((Eigen::MatrixXd(3, 2) << intermediatePositionFrames.front().pos, intermediatePositionFrames.back().pos).finished().transpose(), weights);
    */
    compressedPositionFrames.emplace_back(intermediatePositionFrames.back().time, (Eigen::VectorXd(6) << iFit.row(1).transpose(), intermediatePositionFrames.back().pos).finished());

    intermediatePositionFrames.erase(intermediatePositionFrames.begin(), intermediatePositionFrames.end() - 1);

    // add 2 frames here:
    contactFrames += 2;    
  }
  else if(std::abs(positionError) > pThreshold)
  {
    // upon error bounds tripping, add the f-1 keyframe:
    compressedPositionFrames.push_back(cpPosFrame);
    compressed = 1;

    // remove the irrelevant frames
    intermediatePositionFrames.erase(intermediatePositionFrames.begin(), intermediatePositionFrames.end() - 2);
    
    // add 1 error based frame:
    errorFrames += 1;
  }
  
  // update the previous frames
  cppPosFrame = cpPosFrame;  
  cpPosFrame = cPosFrame;

  if(compressed)
  {
    // reset the error calculations
    positionError = 0;
    positionErrorPrev = 0;
    positionErrorDiffPrev = 0;
    
    qSolve.resetFit(&intermediatePositionFrames);
  }
  return compressed;
}

bool Compressor::compressFrameRotation(bool force)
{  
  // the rotation step:
  // Eigen::VectorXd frame;
  cRotFrame = CompressionEngine::CompressedRotationFrame(intermediateRotationFrames.back().time);
  rotationError = 0.0;
  bool compressed = 0;
  
  if(intermediateRotationFrames.size() >= 2)
    rotationError = RI.linearFit(intermediateRotationFrames, cRotFrame);
  else
    cRotFrame.rot = (Eigen::Vector4d() << 0, intermediateRotationFrames.back().rot.x(), intermediateRotationFrames.back().rot.y(), intermediateRotationFrames.back().rot.z()).finished();

  // keyframe placement:
  if(force || compressedRotationFrames.size() == 0)
  {
    compressedRotationFrames.push_back(cRotFrame);
    compressed = 1;
    
    // remove the frames from the window
    intermediateRotationFrames.erase(intermediateRotationFrames.begin(), intermediateRotationFrames.end() - 1);
  }
  else if(peakDetected(rPThreshold, rotationError, rotationErrorPrev, rotationErrorDiffPrev))
  {
    compressedRotationFrames.push_back(cppRotFrame);
    compressed = 1;

    // remove the frames from the window:
    intermediateRotationFrames.erase(intermediateRotationFrames.begin(), intermediateRotationFrames.end() - 3);

    // Something weird happening if dropping a frame too close. Skip for now
    
    compressedRotationFrames.emplace_back(intermediateRotationFrames.back().time);
    RI.linearFit(intermediateRotationFrames, compressedRotationFrames.back());
    intermediateRotationFrames.erase(intermediateRotationFrames.begin(), intermediateRotationFrames.end() - 1);
    
    contactFlag += 2;
  }
  else if(std::abs(rotationError) > rThreshold)
  {
    compressedRotationFrames.push_back(cpRotFrame);
    compressed = 1;

    // remove the frames from the window
    intermediateRotationFrames.erase(intermediateRotationFrames.begin(), intermediateRotationFrames.end() - 2);
    errorFrames += 1;
  }

  cppRotFrame = cpRotFrame;
  cpRotFrame = cRotFrame;

  if(compressed)
  {
    // reset the errors
    rotationError = 0;
    rotationErrorPrev = 0;
    rotationErrorDiffPrev = 0;
  }

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

  if(std::abs(currentDiff - prevErrorDiff) > threshold)
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
  for(uint32_t i = 0; i < compressedPositionFrames.size(); i++)
  {
    cMat(0, i) = compressedPositionFrames[i].time;
    cMat.col(i).tail(6) = compressedPositionFrames[i].pos;
  }

  return cMat;
}

Eigen::MatrixXd Compressor::getCompressedMatrixRotationRepresentation()
{
  Eigen::MatrixXd cMat(5, compressedRotationFrames.size());
  for(uint32_t i = 0; i < compressedRotationFrames.size(); i++)
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
