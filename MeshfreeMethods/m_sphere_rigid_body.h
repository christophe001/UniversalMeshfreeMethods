/*! \file m_continuous_rigid_body.h */

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

#ifndef _M4_SPHERE_RIGID_BODY_
#define _M4_SPHERE_RIGID_BODY_

#include "m_lagrangian_compute.h"
#include "m_shape.h"

namespace msl {
	class SphereRigidBody : public LagrangianCompute {
	protected:
		double					rho_;
		double					mass_;
		double					r_;
		double					epsilon_;
		friend class ContactMngerSRB;
	public:
		SphereRigidBody(Vec3d center, double r, double rho, Vec3d vel = Vec3d::Zero(), double eps = 0);
		void addAcc(const Vec3d acc) { acc_ += acc; }
		double getMass() const { return mass_; }
		double getRadius() const { return r_; }
		double getRho() const { return rho_; }
		virtual ~SphereRigidBody() {}
		Vec3d	pos_;
		Vec3d	vel_;
		Vec3d	acc_;
		std::string format() const override {
			return "a continuous rigid body!\n";
		}
		void computeForces() override;
		Vec3d getNormal(const Vec3d& pos) { return (pos - pos_).normalized(); }
		double isWithin(const Vec3d& pos) { return (pos - pos_).norm() < r_ + epsilon_; }
	};

}

#endif // !_M4_CONTINUOUS_RIGID_BODY_

