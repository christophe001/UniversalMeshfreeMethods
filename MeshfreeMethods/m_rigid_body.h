/*! \file m_rigid_body.h */

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

#ifndef _M4_RIGID_BODY_
#define _M4_RIGID_BODY_

#include "m_lagrangian_compute.h"

namespace msl {
	class RigidBody : public LagrangianCompute {
	protected:
		Mat3d	moi_;
		Vec3d	center_;
		Vec3d	center_vel_;
		Vec3d	center_acc_;
		Vec3d	omega_;
		Vec3d	omega_dot_;
		void	computeCenter();
		void	computeOmega();
	public:
		RigidBody(std::shared_ptr<SortEnsemble> sorted_ptr);
		virtual ~RigidBody() {}
		std::string format() const override {
			return "a rigid body\n";
		}
		void computeForces() override;
	};
}
#endif