/*
  Position handling for the compression shceme

  Timothy Jeruzalski - 2018
*/
#ifndef QUADRATIC_H
#define QUADRATIC_H

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "CompressionEngine.h"

class Quadratic
{
public:
  Quadratic()
  {
    // set the constraints for the optimization
    constraints = Eigen::MatrixXd(2,3);
    constraints << 0, 0, 1,   1, 1, 1;
    resetFit(NULL);
  }
  
  bool resetFit(std::vector<CompressionEngine::Frame>* intermediateFrames)
  {
    intermediate = intermediateFrames;
    /* FAST methods are only to be used per-object compression */
    fastTTT = Eigen::MatrixXd::Zero(3,3);
    fastTTX = Eigen::MatrixXd::Zero(3,3);
    fastXTX = Eigen::MatrixXd::Zero(3,3);
    NF = Eigen::MatrixXd::Zero(3,3);

    /* Add the previously entered intermediate frames into the running fit alg */
    for(uint32_t j = 0; intermediateFrames != NULL && j < intermediateFrames->size(); j++)
    {
      fastUpdate((*intermediateFrames)[j], (double)j);
    }
    return 1;
  }
  
  bool precompute(std::vector<CompressionEngine::Frame>& intermediateFrames);
  bool fastPrecompute(std::vector<CompressionEngine::Frame>& intermediateFrames);
  
  bool solve(const Eigen::MatrixXd & BCs, Eigen::MatrixXd & weights);
  bool fastSolve(const Eigen::MatrixXd & BCs, Eigen::MatrixXd & weights);

  // do a fast error estimate with partially updated LS fit data
  double fastErrorEstimate(const Eigen::MatrixXd& weights);

private:
  /*
    QUADRATIC SCHEME:
    
     computes the energy in the state from the intermediate frames vector
     TMat is a stacked time matrix of the form:
     [ t₀²  t₀  1 ]
     [ t₁²  t₁  1 ]
     [     ....   ]
     Where the actual time values are normalized such that t₀ = 0 and tₘ = 1 

   */
  double quadraticEnergy(const Eigen::MatrixXd & T, const Eigen::MatrixXd & X, const Eigen::MatrixXd & A);

  // performs minimization in order to solve for the coefficients
  bool determineQuadraticCoefficients(const Eigen::MatrixXd & X, const Eigen::MatrixXd & T, const Eigen::MatrixXd & BCs, Eigen::MatrixXd & A);
  bool fastDetermineQuadraticCoefficients(const Eigen::MatrixXd& TTX, const Eigen::MatrixXd& BCs, Eigen::MatrixXd& A);

  bool fastUpdate(const CompressionEngine::Frame& newestFrame);
  bool fastUpdate(const CompressionEngine::Frame& newestPos, double frameNumber);
  
  // constructs a full matrix from the intermediate frame vector
  bool constructIntermediateStateMatrix(const std::vector<CompressionEngine::Frame>& intermediateFrames, Eigen::MatrixXd& states);

  bool constructTimeMatrix(int32_t steps, Eigen::MatrixXd & t);

  bool constructInvLHS(const Eigen::MatrixXd & T);

  Eigen::MatrixXd invLhs;
  Eigen::MatrixXd constraints; 
  Eigen::MatrixXd timeMatrix;

  // fast error estimate data
  Eigen::MatrixXd fastTTT;
  Eigen::MatrixXd fastTTX;
  Eigen::MatrixXd NF;
  Eigen::MatrixXd fastXTX;
  
  // the intermediate frames
  std::vector<CompressionEngine::Frame>* intermediate = NULL;
};

#endif //QUADRATIC_H
