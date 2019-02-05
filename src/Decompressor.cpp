/*
  Decompressor Code

  Timothy Jeruzalski - 2018
*/

#include "Decompressor.h"

Decompressor::Decompressor()
{
  rData = NULL;
  pData = NULL;
}

bool Decompressor::setCompressedData(Eigen::MatrixXd* rotationData, Eigen::MatrixXd* positionData)
{
  rData = rotationData;
  pData = positionData;

  return 1;
}

CompressionEngine::Frame Decompressor::stateAtTime(double time)
{
  // find the indices and copy out the decompressed data from there
  CompressionEngine::Frame decompTime(time);

  decompTime.rot = decodeRotation(time, findIndexFromTime(rData, time));
  decompTime.pos = decodePosition(time, findIndexFromTime(pData, time));

  return decompTime;
}

Eigen::Quaterniond Decompressor::decodeRotation(double time, int32_t index)
{
  // if before the first time, then return the first rotation
  if(index == 0)
    return completeQuat((*rData).col(0).tail<4>());
  
  // get the normalized time:
  double normT = (time - (*rData)(0, index - 1)) / ((*rData)(0, index) - (*rData)(0, index - 1));

  Eigen::Quaterniond result;
  Eigen::Quaterniond qInit = completeQuat((*rData).col(index - 1).tail<4>());
  Eigen::Quaterniond qFinal = completeQuat((*rData).col(index).tail<4>());

  double alpha = (*rData)(1, index);

  if(std::abs(alpha) < 0.1)
    result = qInit.slerp(normT, qFinal);
  else
  {
    Eigen::AngleAxis<double> diff(qFinal * qInit.conjugate());
    diff.angle() = alpha * normT;

    result = Eigen::Quaterniond(diff) * qInit;
  }

  return result.normalized();
}

Eigen::Vector3d Decompressor::decodePosition(double time, int32_t index)
{
  // if before or is the keyframe, return the first position
  if(index == 0)
    return pData->col(0).segment<3>(3);

  double normT = (time - (*pData)(0, index - 1)) / ((*pData)(0, index) - (*pData)(0, index - 1));

  Eigen::Vector3d result;

  result =
    pData->col(index - 1).segment(3, 3) +
    pData->col(index).segment(0, 3) * normT +
    (pData->col(index).segment(3, 3) - pData->col(index).segment(0, 3) - pData->col(index - 1).segment(3, 3)) * normT * normT;

  return result;
}

Eigen::Quaterniond Decompressor::completeQuat(const Eigen::Vector4d& q)
{
  double wSq = (1 - q(3)*q(3) - q(2)*q(2) - q(1)*q(1));
  return Eigen::Quaterniond(wSq <= 0.001 ? 0 : std::sqrt(wSq), q(1), q(2), q(3)).normalized();
}

void Decompressor::flipQuats(Eigen::Quaterniond& q1, Eigen::Quaterniond& q2)
{
  if(q1.coeffs().dot(q2.coeffs()) < 0)
    q2.coeffs() *= -1;

  // just in case
  q1.normalize();
  q2.normalize();
}

int32_t Decompressor::findIndexFromTime(const Eigen::MatrixXd* data, double time)
{
  if((*data).cols() < 2)
  {
    std::cout << "Too few keyframes" << std::endl;
    return -1;
  }

  int32_t l_index = 0, u_index = data->cols() - 1;

  if(time < (*data)(0, 0))
    return l_index;
  else if(time > (*data)(0, data->cols() - 1))
    return u_index;

  while((u_index - l_index) > 1)
  {
    int32_t m_index = l_index + std::floor((u_index - l_index) / 2.0);

    if((*data)(0, m_index) < time)
      l_index = m_index;
    else
      u_index = m_index;
  }

  return u_index;
}
