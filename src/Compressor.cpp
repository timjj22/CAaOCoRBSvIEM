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

bool Compressor::compressFrame(Frame& newFrame)
{
  windowFrames.push_back(newFrame);

  // the position step:
  qSolve.fastPrecompute(intermediateFrames);
  Eigen::MatrixXd weights;
  qSolve.fastSolve((Eigen::MatrixXd(3, 2) << intermediateFrames.front().pos, intermediateFrames.back().pos).finished().transpose(), weights);
  positionError = qSolve.fastErrorEstimate(weights);
  // now set the current compressed frame
  cPosFrame = CompressionEngine::CompressedPositionFrame(newFrame.time, (Eigen::VectorXd(6) << weights.row(1).transpose(), intermediateFrames.back().pos).finished());

  // the rotation step:
  Eigen::VectorXd frame;
  if(intermediateFrames.size() >= 2)
    rotationError = RI.linearFit(intermediateFrames, frame);
  else
    frame = (Eigen::VectorXd(4) << 0, newFrame.rot.x(), newFrame.rot.y(), newFrame.rot.z()).finished();

  // check the error conditions
  double errorDiff = positionError - positionErrorPrev;
  if(std::abs(errorDiff - positionErrorDiffPrev) > pPThreshold)
  positionErrorDiffPrev = errorDiff;

  // update the previous frames
  if(cpPosFrame.size() > 0)
    cppPosFrame = cpPosFrame;
  if(cpRotFrame.size() > 0)
    cppRotFrame = cpRotFrame;

  // swap out the old frames
  if(pFrame.size() > 0)
    ppFrame = pFrame;
  if(newFrame.size() > 0)
    pFrame = newFrame;
}

Eigen::MatrixXd getCompressedMatrixPositionRepresentation()
{
  Eigen::MatrixXd cMat(7, compressedPositionFrames.size());
  for(int32_t i = 0; i < CompressedPositionFrames.size(); i++)
  {
    cMat(0, i) = compressedPositionFrames[i].time;
    cMat.col(i).tail(6) = compressedPositionFrames[i].pos;
  }

  return cMat;
}

Eigen::MatrixXd getCompressedMatrixRotationRepresentation()
{
  Eigen::MatrixXd cMat(5, compressedRotationFrames.size());
  for(int32_t i = 0; i < CompressedRotationFrames.size(); i++)
  {
    cMat(0, i) = compressedRotationFrames[i].time;
    cMat.col(i).tail(4) = compressedRotationFrames[i].rot;
  }

  return cMat;
}

std::vector<CompressionEngine::CompressedPositionFrame> getCompressedVectorPositionRepresentation()
{
  return compressedPositionFrames;
}

std::vector<CompressionEngine::CompressedRotationFrame> getCompressedVectorRotationRepresentation()
{
  return compressedRotationFrames;
}
