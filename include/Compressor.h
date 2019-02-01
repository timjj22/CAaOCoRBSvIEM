/*
  Compressor code for the 2018 SCA paper "Collision-Aware and Online Compression of Rigid Body Simulations via Integrated Error Minimization"

  Timothy Jeruzalski - 2018
 */
#ifndef COMPRESSOR_H
#define COMPRESSOR_H

#include "CompressionEngine.h"

class Compressor : CompressionEngine {
public:  
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Compressor(double positionThreshold = 5e-2,
	     double positionPeakThreshold = 5e-4,
	     double rotationThreshold = 1e-2,
	     double rotationPeakThreshold = 5e-3);
  ~Compressor();

  // add a frame to the compressed data
  bool compressFrame(const Frame& newFrame, bool force = false);

  bool finalizeRepresentation(const Frame& newFrame);
  bool finalizeRepresentation();

  Eigen::MatrixXd getCompressedMatrixPositionRepresentation();
  Eigen::MatrixXd getCompressedMatrixRotationRepresentation();
  std::vector<CompressionEngine::CompressedPositionFrame> getCompressedVectorPositionRepresentation();
  std::vector<CompressionEngine::CompressedRotationFrame> getCompressedVectorRotationRepresentation();
  
private:
  // using the intermediateFrames, compute the error and compressed frames
  bool compressFramePosition();
  bool compressFrameRotation();

  bool peakDetected(const double& error, double& error);
  
  // the fast position fit and error metric solver
  Quadratic qSolve;
  
  double pThreshold;   // Position
  double pPThreshold;  // Position Peak
  double rThreshold;   // Rotation  
  double rPThreshold;  // Rotation Peak
  
  // peak threshold calculations
  double positionError;
  double positionErrorPrev;
  double positionErrorDiffPrev;

  double rotationError;
  double rotationErrorPrev;
  double rotationErrorDiffPrev;

  // previous and previous-previous compressed frames for use with the peak detection phase
  CompressionEngine::CompressedPositionFrame cPosFrame;
  CompressionEngine::CompressedPositionFrame cpPosFrame;
  CompressionEngine::CompressedPositionFrame cppPosFrame;
  
  CompressionEngine::CompressedRotationFrame cRotFrame;
  CompressionEngine::CompressedRotationFrame cpRotFrame;
  CompressionEngine::CompressedRotationFrame cppRotFrame;

  // previous and previous-previous frame
  CompressionEngine::Frame pFrame;
  CompressionEngine::Frame ppFrame;
  
  // The whole amount of sliding window frames does not need to be kept. Only a small subset of frames are needed for the rotation fitting
  std::vector<CompressionEngine::Frame> intermediatePositionFrames;
  std::vector<CompressionEngine::Frame> intermediateRotationFrames;
  
  // compressed states
  std::vector<CompressionEngine::CompressedPositionFrame> compressedPositionFrames;
  std::vector<CompressionEngine::CompressedRotationFrame> compressedRotationFrames;
  
  // stats gathering
  int32_t contactFrames = 0;
  int32_t errorFrames = 0;
  bool contactFlag = 0;
}

#endif//COMPRESSOR_H
