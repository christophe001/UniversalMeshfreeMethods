/*! \file m_contact_manager.h */

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
#ifndef _M4_CONTACT_MANAGER_
#define _M4_CONTACT_MANAGER_

#include "m_types.h"
#include "m_ensemble.h"
#include "m_sort_ensemble.h"

namespace msl {
	class ContactManager : public MsObj {
	protected:
		std::shared_ptr<SortEnsemble>	slave_;
		std::shared_ptr<SortEnsemble>	master_;
		std::shared_ptr<Ensemble>		ensemble_s_;
		std::shared_ptr<Ensemble>		ensemble_m_;
		double							epsilon_;
		Vec3d							contact_pmin_;
		Vec3d							contact_pmax_;
		Vec3i							cell_min_;
		Vec3i							cell_max_;
		Vec3d*							spos_;
		Vec3d*							mpos_;
		Vec3d*							sacc_;
		Vec3d*							macc_;
		DomainConfig*					domain_cfg_;
		double							dv_;
		double							force_const_;
		double							time_to_contact_;

	public:
		ContactManager(std::shared_ptr<SortEnsemble> master,
			std::shared_ptr<SortEnsemble> slave, double epsilon);
		void setEpsilon(const double& epsilon) { epsilon_ = epsilon; }
		void setForceParams(const double& dv, const double& fc) { dv_ = dv; force_const_ = fc; }
		ContactManager(const ContactManager& cm) = delete;
		ContactManager& operator=(const ContactManager& cm) = delete;
		virtual ~ContactManager() {}
		typedef void (ContactManager::*handler)(int i, int j);
		void computeContact(handler hptr = &penaltyHandler);
		void computeContactTime();
		void updateContactZone();
		void penaltyHandler(int i, int j);
		double getTimeToContact() { return time_to_contact_; }
	};
}
#endif // !_M4_CONTACT_MANAGER_
