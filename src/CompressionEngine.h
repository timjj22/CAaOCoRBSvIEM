/*
  CompressionEngine for shared functionality with compressing and decompressing 

  Timothy Jeruzalski - 2018
 */
#ifndef COMPRESSIONENGINE_H
#define COMPRESSIONENGINE_H

#include <Eigen/Core.h>

class CompressionEngine {
public:
  struct Frame {
    Frame(double t, Eigen::Quaterniond rotation, Eigen::Vector3d position) : time(t), rot(rotation), pos(position) { };
    ~Frame() { };

    double time;
    Eigen::Quaterniond rot;
    Eigen::Vector3d pos;
  };
  struct CompressedPositionFrame {
    CompressedFrame(double t, Eigen::VectorXd position) : time(t), pos(position) { };
    ~CompressedFrame() { };

    double time;
    Eigen::VectorXd pos;
  };
  struct CompressedRotationFrame {
    CompressedFrame(double t, Eigen::Vector4d rotation) : time(t), rot(rotation) { };
    ~CompressedFrame() { };
    
    double time;
    Eigen::Vector4d rot; // is not a quaternion, but 3 quaternion components and an angular velocity
  };
};

#endif//COMPRESSIONENGINE_H
