#include "Quadratic.h"


double Quadratic::quadraticEnergy(const Eigen::MatrixXd & T, const Eigen::MatrixXd & X, const Eigen::MatrixXd & A)
{
  /*
    The goal of the optimization is to solve the following where x(t) ϵ ℝ³ is the data term and α, β, γ ϵ ℝ³
    are the unknowns to find

    min_{α, β, γ} ½ ∑‖ x(t) - (αt² + βt + γ) ‖² ∆t

    expanding we can for an energy in terms of A = [α, β, γ]ᵀ, T is stacked time to various powers, X is stacked states

    min_{α, β, γ} ½ ‖ X - TA ‖² ∆t

    For which we can expand the Frobenius norm and express an energy:
    E(A) = ½tr(XᵀX) - tr(XᵀTA) + ½tr(AᵀTᵀTA)

   */

  return 0.5*(X.transpose()*X).trace() - (X.transpose()*T*A).trace() + 0.5*(A.transpose()*T.transpose()*T*A).trace();

}

bool Quadratic::constructInvLHS(const Eigen::MatrixXd & TTT)
{
  /*    
    adding constraints to the optimization through lagrange multipliers

    [ TᵀT   Cᵀ ] [ A ]   [    XᵀT   ]
    [ C     0  ] [ λ ] = [ [x₀ xₘ]ᵀ ]

    which can be solved by inverting the LHS, and can be applied for all objects in the scene for current timestep
    
    (TᵀT is TTT):
   */

  Eigen::MatrixXd LHS = Eigen::MatrixXd::Zero(TTT.cols() + constraints.rows(), TTT.cols() + constraints.rows());

  LHS.block(0, 0, TTT.cols(), TTT.cols()) = TTT;
  LHS.block(TTT.cols(), 0, constraints.rows(), constraints.cols()) = constraints;
  LHS.block(0, TTT.cols(), constraints.cols(), constraints.rows()) = constraints.transpose();

  invLhs = LHS.inverse();

  // std::cout << invLhs << std::endl << std::endl;
  
  return 1;
}

bool Quadratic::constructTimeMatrix(int32_t steps, Eigen::MatrixXd & t)
{
  /*
    the time matrix for quadratic energy computations is of the form:
    [ t₀²  t₀  1 ]
    [ t₁²  t₁  1 ]
    [     ....   ]
    
    but is normalized between 0 and 1
  */

  t = Eigen::MatrixXd(steps, 3);
  for(int32_t i = 0; i < steps; i++)
  {
    double time = steps > 1 ? i /((double)steps - 1) : 0;

    t(i, 0) = time*time;
    t(i, 1) = time;
    t(i, 2) = 1;
  }

  return 1;
}

bool Quadratic::determineQuadraticCoefficients(const Eigen::MatrixXd & X, const Eigen::MatrixXd & T, const Eigen::MatrixXd & BCs, Eigen::MatrixXd & A)
{
  return fastDetermineQuadraticCoefficients(T.transpose()*X, BCs, A);
}

bool Quadratic::fastDetermineQuadraticCoefficients(const Eigen::MatrixXd& TTX, const Eigen::MatrixXd& BCs, Eigen::MatrixXd& A)
{
  // determine the coefficients with the given TᵀT and TᵀX
  // will determine the optimal set of parameters for the fitting for a single object
  Eigen::MatrixXd rhs(TTX.cols() + BCs.rows(),TTX.cols());  
  
  rhs.block(0,0, TTX.cols(), TTX.cols()) = TTX;
  rhs.block(TTX.cols(), 0, BCs.rows(), BCs.cols()) = BCs;

  Eigen::MatrixXd temp = invLhs*rhs;

  //std::cout << rhs  << std::endl << std::endl;
  //std::cout << temp << std::endl << std::endl;

  A = temp.block(0, 0, TTX.cols(), TTX.cols());
  return 1;      
}

bool Quadratic::constructIntermediateStateMatrix(const std::vector<CompressionEngine::Frame>& intermediateFrames, Eigen::MatrixXd & states)
{
  // cannot make if no data
  if(intermediateFrames.size() == 0)
    return 0;

  // can only deal with single object at a time
  states = Eigen::MatrixXd(intermediateFrames.size(), 3);
  
  bool success = 1;
  // fill with just the position data from each frame
  for(uint32_t i = 0; i < intermediateFrames.size(); i++)
  {
    (states.block(i, 0, 1, 3)) = intermediateFrames[i].pos.transpose();
  }
  return success;
}

bool Quadratic::precompute(std::vector<CompressionEngine::Frame>& intermediateFrames)
{  
  if(!constructTimeMatrix(intermediateFrames.size(), timeMatrix))
    printf("ERROR creating quadratic time matrix\n");
  
  if(timeMatrix.rows() == 0)
    return 1;
  
  if(!constructInvLHS(timeMatrix.transpose()*timeMatrix))
    printf("ERROR creating inverted quadratic matrix\n");

  intermediate = &intermediateFrames;
  return 1;
}

bool Quadratic::fastPrecompute(std::vector<CompressionEngine::Frame>& intermediateFrames)
{
  // TᵀT can be precomputed for all of the objects in the scene
  intermediate = &intermediateFrames;
  
  if(intermediate->size() < 1)
    return 0;

  double nf = (double)intermediate->size() - 1;
  
  NF <<
    (1/pow(nf, 2)), 0, 0,
    0, (1/nf), 0,
    0, 0, 1;

  fastTTT <<
    (1/30.)*(1/pow(nf, 4))*nf*(nf + 1)*(2*nf + 1)*(3*pow(nf, 2) + 3*nf - 1),   (1/4.)*(1/pow(nf, 3))*pow(nf, 2)*pow(nf + 1, 2),    (1/pow(nf, 2))*(1/6.)*(nf)*(nf+1)*(2*nf+1),
    (1/4.)*(1/pow(nf, 3))*pow(nf, 2)*pow(nf + 1, 2),                           (1/pow(nf, 2))*(1/6.)*(nf)*(nf+1)*(2*nf+1),         (1/nf)*(1/2.)*(nf)*(nf + 1),
    (1/pow(nf, 2))*(1/6.)*(nf)*(nf+1)*(2*nf+1),                                (1/nf)*(1/2.)*(nf)*(nf + 1),                        nf + 1.0;

  // Eigen::MatrixXd TtT;
  // constructTimeMatrix(intermediateFrames.size(), TtT);
  // constructInvLHS(TtT);
  // std::cout << "SLOW: \n" << TtT.transpose()*TtT << std::endl;
  
  constructInvLHS(fastTTT);
  
  return 1;
}

bool Quadratic::fastUpdate(const CompressionEngine::Frame& newestFrame)
{
  return fastUpdate(newestFrame, (double)intermediate->size() - 1.0);
}

bool Quadratic::fastUpdate(const CompressionEngine::Frame& newestFrame, double frameNumber)
{
  /*
    This will update the data dependent terms in the fast update least squares solve
  */  
  fastTTX = fastTTX + (Eigen::VectorXd(3) << (pow(frameNumber, 2)), (frameNumber), 1).finished()*newestFrame.pos.transpose();
  fastXTX = fastXTX + newestFrame.pos*newestFrame.pos.transpose();
  return 1;
}

bool Quadratic::fastSolve(const Eigen::MatrixXd & BCs, Eigen::MatrixXd & weights)
{
  if(intermediate->size() < 1)
    return 0;

  CompressionEngine::Frame newestPos = intermediate->back();
  fastUpdate(newestPos);

  // do a linear fit if there aren't enough frames to do a quadratic fit
  if(intermediate->size() < 3)
  {
    weights = (Eigen::MatrixXd(3,3) << 0,0,0, BCs.row(1).x() - BCs.row(0).x(), BCs.row(1).y() - BCs.row(0).y(), BCs.row(1).z() - BCs.row(0).z(), BCs.row(0).x(), BCs.row(0).y(), BCs.row(0).z() ).finished();
    //std::cout << weights << std::endl;
  }
  else
    fastDetermineQuadraticCoefficients(NF*fastTTX, BCs, weights);
  // std::cout << weights << std::endl << std::endl;
  return 1;
}

double Quadratic::fastErrorEstimate(const Eigen::MatrixXd& weights)
{
  if(intermediate->size() <= 3)
    return 0;
  else
    return (fastXTX +
	    weights.transpose()*fastTTT*weights -
	    2*weights.transpose()*NF*fastTTX).trace();
}

bool Quadratic::solve(const Eigen::MatrixXd & BCs, Eigen::MatrixXd & weights)
{
  // // update the data for the fast error estimate
  // Eigen::VectorXd newestPos;
  // // get the number of frames in the current window (used for the fast calculation)
  // double nf = (double)intermediate->size() - 1;

  // if(intermediate->size() >= 1)
  // {
  //   NF << (1/pow(nf, 2)), 0, 0,
  //     0, (1/nf), 0,
  //     0, 0, 1;
  //   dataFormat.getPositionFromState(&(intermediate->back()), object, newestPos);
  //   // due to the normalized time, these values are determined through series definitions
  //   fastTTT <<
  //     (1/30.)*(1/pow(nf, 4))*nf*(nf + 1)*(2*nf + 1)*(3*pow(nf, 2) + 3*nf - 1),   (1/4.)*(1/pow(nf, 3))*pow(nf, 2)*pow(nf + 1, 2),    (1/pow(nf, 2))*(1/6.)*(nf)*(nf+1)*(2*nf+1),
  //     (1/4.)*(1/pow(nf, 3))*pow(nf, 2)*pow(nf + 1, 2),                           (1/pow(nf, 2))*(1/6.)*(nf)*(nf+1)*(2*nf+1),         (1/nf)*(1/2.)*(nf)*(nf + 1),
  //     (1/pow(nf, 2))*(1/6.)*(nf)*(nf+1)*(2*nf+1),                                (1/nf)*(1/2.)*(nf)*(nf + 1),                        nf + 1.0;
  //   fastTTX = fastTTX + (Eigen::VectorXd(3) << (pow(nf, 2)), (nf), 1).finished()*newestPos.transpose();
  //   fastXTX = fastXTX + newestPos*newestPos.transpose();
  // }
  // if not enough time history, then do a linear fit (otherwise gets a bunch of NaN's)
  if(timeMatrix.rows() < 3)
  {
    weights = (Eigen::MatrixXd(3,3) << 0,0,0, BCs.row(1).x() - BCs.row(0).x(), BCs.row(1).y() - BCs.row(0).y(), BCs.row(1).z() - BCs.row(0).z(), BCs.row(0).x(), BCs.row(0).y(), BCs.row(0).z() ).finished();
    //std::cout << weights << std::endl;
    return 1;
  }

  // after precomputing, the solve is the same between objects (with different RHS) 
  Eigen::MatrixXd states;
  constructIntermediateStateMatrix(*intermediate, states);

  determineQuadraticCoefficients(states, timeMatrix, BCs, weights);
  //fastXTX = states.transpose()*states;
  return 1;
}
