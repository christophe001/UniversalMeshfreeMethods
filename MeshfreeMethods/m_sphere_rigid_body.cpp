/*! \file m_continuous_rigid_body.cpp */

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

#include "m_sphere_rigid_body.h"

namespace msl {
	
	SphereRigidBody::SphereRigidBody(Vec3d center, double r, double rho, Vec3d vel, double eps) {
		pos_ = center;
		r_ = r;
		rho_ = rho;
		epsilon_ = eps;
		vel_ = vel;
		acc_ = Vec3d::Zero();
		mass_ = 4.0 / 3.0 * pi * r_ * r_ * r_ * rho_;
		class_name_ = "SphereRigidBody";
	}


	void SphereRigidBody::computeForces() {
	
	}

}