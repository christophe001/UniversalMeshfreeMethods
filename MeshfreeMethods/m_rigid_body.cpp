/*! \file m_rigid_body.cpp */

//@HEADER
// ************************************************************************
//
//                         M4_Technologies
//                 Copyright (2017) by Wentao Xu
//
// M4 stands for Multiscale Mesh-based Meshfree Method
// Any questions, please contact:
// Wentao Xu   wx2151@columbia.edu
//
// ************************************************************************
//@HEADER

#include "m_rigid_body.h"

namespace msl {
	void RigidBody::computeCenter() {
		center_ = Vec3d::Zero();
		center_vel_ = Vec3d::Zero();
		center_acc_ = Vec3d::Zero();
		for (int i = 0; i < np_; i++) {
			//! center of mass
			center_ += pos_[i] / np_;
			//! center velocity by conservation of momemtum
			center_vel_ += vel_[i] / np_;
			//! center acceleration by conservation of force
			center_acc_ += acc_[i] / np_;
		}
	}

	void RigidBody::computeOmega() {
		Vec3d torque = Vec3d::Zero();
		Vec3d am = Vec3d::Zero();
		double mass = ensemble_ptr_->getMass();
		for (int i = 0; i < np_; i++) {
			torque += (pos_[i] - center_).cross(acc_[i]) * mass;
			am += (pos_[i] - center_).cross(vel_[i] - center_vel_) * mass;
		}
		omega_dot_ = moi_.inverse() * torque;
		omega_ = moi_.inverse() * am;
	}

	RigidBody::RigidBody(std::shared_ptr<SortEnsemble> sorted_ptr)
		: LagrangianCompute(sorted_ptr) {
		class_name_ = "RigidBody";
		double mass = ensemble_ptr_->getMass();
		computeCenter();
		moi_ = Mat3d::Zero();
		for (int i = 0; i < np_; i++) {
			Vec3d vec = pos_[i] - center_;
			moi_(0, 0) += mass * (vec(1) * vec(1) + vec(2) * vec(2));
			moi_(1, 1) += mass * (vec(0) * vec(0) + vec(2) * vec(2));
			moi_(2, 2) += mass * (vec(0) * vec(0) + vec(1) * vec(1));
			moi_(0, 1) -= mass * vec(0) * vec(1);
			moi_(1, 0) -= mass * vec(0) * vec(1);
			moi_(0, 2) -= mass * vec(0) * vec(2);
			moi_(2, 0) -= mass * vec(0) * vec(2);
			moi_(1, 2) -= mass * vec(1) * vec(2);
			moi_(2, 1) -= mass * vec(1) * vec(2);
		}
	}

	void RigidBody::computeForces() {
		computeCenter();
		computeOmega();
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++)
			acc_[i] = center_acc_ + omega_dot_.cross(pos_[i] - center_) 
			+ omega_.cross(vel_[i] - center_vel_);
	}

	void RigidBody::applyForce(const Vec3d& acc) {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++) {
			acc_[i] = acc;
		}
		center_acc_ = acc;
	}

	void RigidBody::applySpeed(const Vec3d & vel) {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++)
			vel_[i] = vel;
		center_vel_ = vel;
	}
}
