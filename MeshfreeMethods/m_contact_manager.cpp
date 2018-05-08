/*! \file m_contact_manager.cpp */

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

#include "m_contact_manager.h"
#include <omp.h>

namespace msl {

	//==============================================================================
	/// Ctor
	//==============================================================================
	ContactManager::ContactManager(std::shared_ptr<SortEnsemble> master, 
		std::shared_ptr<SortEnsemble> slave, double epsilon) 
		: slave_(slave), master_(master), 
		epsilon_(epsilon), time_to_contact_(0)
	{
#ifdef _DEBUG_
		if (*(slave_->domain_cfg_) != *(master_->domain_cfg_))
			throwException("ContactManager", "Ensemble have different domain configuration");
#endif // _DEBUG_
		domain_cfg_ = master_->getDomainConfig();
		ensemble_s_ = slave_->ensemble_ptr_;
		ensemble_m_ = master_->ensemble_ptr_;
		spos_ = ensemble_s_->getPos();
		mpos_ = ensemble_m_->getPos();
		sacc_ = ensemble_s_->getAcc();
		macc_ = ensemble_m_->getAcc();
		class_name_ = "ContactManager";
		updateContactZone();
		computeContactTime();
	}

	//==============================================================================
	/// Contact compute
	//==============================================================================
	void ContactManager::computeContact(handler hptr) {
		if ((contact_pmax_.array() >= contact_pmin_.array()).all()) {
			int x = cell_max_[0] - cell_min_[0] + 1;
			int y = cell_max_[1] - cell_min_[1] + 1;
			int z = cell_max_[2] - cell_min_[2] + 1;
			int cell_num = x * y * z;
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
			for (int i = 0; i < cell_num; i++) {
				Vec3i cur(i % x, (i % (x * y)) / x, i / (x * y));
				Vec3i cell = cur + cell_min_;
				if (master_->hasCell(cell)) {
					auto it_m = master_->findCellStart(cell);
					auto vectors = slave_->getAllAdjacentCells(cell);
					for (int i = master_->begin(it_m); i < master_->end(it_m); i++) {
						if (!vectors.empty())
						for (auto sc_it : vectors) {
							for (int j = slave_->begin(sc_it); j < slave_->end(sc_it); j++) {
								(this->*hptr)(i, j);
							}
						}
					}
					if (slave_->hasCell(cell)) {
						auto it_s = slave_->findCellStart(cell);
						for (int i = master_->begin(it_m); i < master_->end(it_m); i++) {
							for (int j = slave_->begin(it_s); j < slave_->end(it_s); j++) {
								(this->*hptr)(i, j);
							}
						}
					}
				}					
			}
		}
	}

	//==============================================================================
	/// compute time to contact
	//==============================================================================
	void ContactManager::computeContactTime() {
		Vec2d sv = ensemble_s_->getVelMinMax();
		Vec2d mv = ensemble_m_->getVelMinMax();
		time_to_contact_ =  epsilon_ / (sv[1] + mv[1]);
	}

	//==============================================================================
	/// Update contact zone
	//==============================================================================
	void ContactManager::updateContactZone() {
		Vec3d s_pmin = ensemble_s_->getPosMin();
		Vec3d s_pmax = ensemble_s_->getPosMax();
		Vec3d m_pmin = ensemble_m_->getPosMin();
		Vec3d m_pmax = ensemble_m_->getPosMax();
		Vec3i cell_offset = domain_cfg_->getCellOffset();
		Vec3d aug_m_pmin = m_pmin - epsilon_ * cell_offset.cast<double>();
		Vec3d aug_m_pmax = m_pmax + epsilon_ * cell_offset.cast<double>();
		contact_pmin_ = m_pmin - epsilon_ * cell_offset.cast<double>();
		contact_pmax_ = m_pmax + epsilon_ * cell_offset.cast<double>();
		contact_pmin_ = aug_m_pmin.cwiseMax(s_pmin);
		contact_pmax_ = aug_m_pmax.cwiseMin(s_pmax);
		cell_min_ = domain_cfg_->getCellNum(contact_pmin_);
		cell_max_ = domain_cfg_->getCellNum(contact_pmax_);
	}

	//==============================================================================
	/// Penalty handler
	//==============================================================================
	void ContactManager::penaltyHandler(int i, int j) {
		if ((spos_[j] - mpos_[i]).norm() < epsilon_) {
			Vec3d eta = spos_[j] - mpos_[i];
			Vec3d force = force_const_ * dv_ * eta / eta.norm() * log(epsilon_ / eta.norm());
			macc_[i] -= force;
			sacc_[j] += force;
		}
	}

}