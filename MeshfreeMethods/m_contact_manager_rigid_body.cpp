#include "m_contact_manager_rigid_body.h"
#include <cmath>

namespace msl {
	ContactMngerSRB::ContactMngerSRB(std::shared_ptr<RigidBody>& master,
		std::shared_ptr<SortEnsemble>& slave, double radius, double dp, double force_const)
	{
		class_name_ = "ContactMngerSRB";
		slave_ = slave;
		master_ = master;
		dp_ = dp;
		radius_ = radius;
		force_const_ = force_const;
		np_ = slave_->getNp();
		s_mass_ = slave_->getEnsemble()->getMass();
		m_mass_ = master_->getEnsemble()->getMass() * master_->getEnsemble()->getNp();
		s_pos_ = slave_->getEnsemble()->getPos();
		s_vel_ = slave_->getEnsemble()->getVel();
		s_acc_ = slave_->getEnsemble()->getAcc();
		eps_ = radius_ + dp_ / 2.0;
		if (slave_->getEnsemble()->hasScalarAttribute("damage"))
			damage_ = slave_->getEnsemble()->getScalarAttrPtr("damage")->getAttr();
		else
			damage_ = nullptr;
	}

	void ContactMngerSRB::computeForces() {
		master_->computeCenter();
		center_ = master_->getCenterPos();
		vel_ = master_->getCenterVel();
		acc_ = Vec3d::Zero();
//#ifdef _WITH_OMP_
//#pragma omp parallel for schedule(static)
//#endif // _WITH_OMP
		for (int i = 0; i < np_; i++) {
			Vec3d eta = center_ - s_pos_[i];
			double nm = eta.norm();
			double f = force(nm);
			s_acc_[i] -= f * eta / (s_mass_ * nm);
			acc_ += f * eta / (m_mass_ * nm);
		}
		master_->applyForce(acc_);
	}

	double ContactMngerSRB::force(const double & d) const {
		return d > eps_ ? 0.0 :  d == 0.0 ? force_const_ * 10000.0 : std::min(force_const_ * log(eps_ / d), force_const_ * 10000.0);
	}

	void ContactMngerSRB::enforceNoSlip() {
		if (damage_ == nullptr)
			throwException("enforceNoSlip", "No damage found!");
		master_->computeCenter();
		center_ = master_->getCenterPos();
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
		for (int i = 0; i < np_; i++) {
			Vec3d eta = center_ - s_pos_[i];
			double nm = eta.norm();
			if (nm < eps_ && damage_[i] >= 0.99) {
				s_vel_[i] = master_->getCenterVel();
				s_acc_[i] = master_->getCenterAcc();
			}
		}

	}

	void ContactMngerSRB::rigidBodyEnforceAccVel(const Vec3d & acc, const Vec3d & vel) {
		master_->applyForce(acc);
		master_->applySpeed(vel);
	}

	//void ContactMngerSRB::computeContact(ContactManager::handler handler) {
		
	//}
	
	//void ContactMngerSRB::updateContactZone() {

	//}
}