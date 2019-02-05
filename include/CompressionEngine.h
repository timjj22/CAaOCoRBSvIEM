/*
  CompressionEngine for shared functionality with compressing and decompressing 

  Timothy Jeruzalski - 2018
 */
#ifndef COMPRESSIONENGINE_H
#define COMPRESSIONENGINE_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <iostream>

class CompressionEngine {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  struct Frame {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Frame() : time(0) { rot = Eigen::Quaterniond::Identity(); pos = Eigen::Vector3d::Zero(); };
    Frame(double t) : time(t) { rot = Eigen::Quaterniond::Identity(); pos = Eigen::Vector3d::Zero(); };
    Frame(double t, Eigen::Quaterniond rotation, Eigen::Vector3d position) : time(t), rot(rotation), pos(position) { };
    Frame(double t, Eigen::VectorXd rotation, Eigen::VectorXd position) : time(t),
									  rot(rotation(3), rotation(0), rotation(1), rotation(2)),
									  pos(position(0), position(1), position(2)) {};
    ~Frame() { };

    double time;
    Eigen::Quaterniond rot;
    Eigen::Vector3d pos;
  };
  struct CompressedPositionFrame {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    CompressedPositionFrame() : time(0) { pos = Eigen::VectorXd::Zero(6);}
    CompressedPositionFrame(double t) : time(t) { pos = Eigen::VectorXd::Zero(6); };
    CompressedPositionFrame(double t, Eigen::VectorXd position) : time(t), pos(position) { };
    ~CompressedPositionFrame() { };

    double time;
    Eigen::VectorXd pos;
  };
  struct CompressedRotationFrame {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    CompressedRotationFrame() : time(0) { rot = Eigen::Vector4d::Zero(); };
    CompressedRotationFrame(double t) : time(t) { rot = Eigen::Vector4d::Zero(); };
    CompressedRotationFrame(double t, Eigen::Vector4d rotation) : time(t), rot(rotation) { };
    ~CompressedRotationFrame() { };    
    
    double time;
    Eigen::Vector4d rot; // is not a quaternion, but 3 quaternion components and an angular velocity
  };
};

#endif//COMPRESSIONENGINE_H
