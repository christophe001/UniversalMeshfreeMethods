/*! \file m_contact_manager_rigid_body.h */

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

#ifndef _M4_CONTACT_MANAGER_RIGID_BODY_
#define _M4_CONTACT_MANAGER_RIGID_BODY_

#include "m_contact_manager.h"
#include "m_meshfree_base.h"

namespace msl {

	class ContactMngerSRB : public MsObj {
	protected:
		Vec3d center_;
		Vec3d vel_;
		Vec3d acc_;
		double radius_;
		double dp_;
		double m_mass_;
		double s_mass_;
		double force_const_;
		double eps_;
		int np_;
		double* damage_;
		Vec3d* s_pos_;
		Vec3d* s_acc_;
		Vec3d* s_vel_;
		std::shared_ptr<SortEnsemble> slave_;
		std::shared_ptr<RigidBody> master_;
	public:
		ContactMngerSRB(std::shared_ptr<RigidBody>& master,
			std::shared_ptr<SortEnsemble>& slave, double radius, double dp, double force_const);
		void update() { 
			master_->computeCenter(); 
			center_ = master_->getCenterPos(); 
			vel_ = master_->getCenterVel();
			acc_ = master_->getCenterAcc();
		}
		Vec3d getCenterPos() const { return center_; }
		Vec3d getCenterVel() const { return vel_; }
		Vec3d getCenterAcc() const { return acc_; }
		double getRadius() const { return radius_; }
		double getDp() const { return dp_; }
		double getMasterTotalMass() const { return m_mass_; }
		double getSlaveSingleMass() const { return s_mass_; }
		double getEpsilon() const { return eps_; }
		~ContactMngerSRB() {}
		//void computeContact(ContactManager::handler handler 
			//= static_cast<ContactManager::handler>(&doNothing)) override;
		//void updateContactZone() override;
		void computeForces();
		//void updateContactZone() override;
		double force(const double& nm) const;

		void enforceNoSlip();
		void rigidBodyEnforceAccVel(const Vec3d& acc, const Vec3d& vel);
	};

}


#endif // !_CONTACT_MANAGER_RIGID_BODY_

