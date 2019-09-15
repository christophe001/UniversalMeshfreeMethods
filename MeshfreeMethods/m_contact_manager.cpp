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
#include <algorithm>
#include <iostream>

namespace msl {

	//==============================================================================
	/// Ctor
	//==============================================================================
	ContactManager::ContactManager(std::shared_ptr<SortEnsemble> master, 
		std::shared_ptr<SortEnsemble> slave, double epsilon, double dp,  double dt) 
		: slave_(slave), master_(master)
	{
#ifdef _DEBUG_
		if (*(slave_->domain_cfg_) != *(master_->domain_cfg_))
			throwException("ContactManager", "Ensemble have different domain configuration");
#endif // _DEBUG_
		setEpsilonDtMax(epsilon, dt);
		dp_ = dp;
		ensemble_s_ = slave_->ensemble_ptr_;
		ensemble_m_ = master_->ensemble_ptr_;
		spos_ = ensemble_s_->getPos();
		mpos_ = ensemble_m_->getPos();
		svel_ = ensemble_s_->getVel();
		mvel_ = ensemble_m_->getVel();
		sacc_ = ensemble_s_->getAcc();
		macc_ = ensemble_m_->getAcc();
		domain_cfg_ = DomainConfig(*(master_->getDomainConfig()));
		contact_pmin_ = Vec3d::Ones();
		contact_pmax_ = -1.0 * Vec3d::Ones();
		class_name_ = "ContactManager";
		updateContactZone();
	}

	void ContactManager::setEpsilonDtMax(const double & epsilon, const double & dt_max) {
		epsilon_ = epsilon;
		dt_ = dt_max;
		dt_max_ = dt_max;
	}

	//==============================================================================
	/// Contact compute
	//==============================================================================
	void ContactManager::computeContact(handler hptr) {
		if ((contact_pmax_.array() > contact_pmin_.array()).all()) {
			int x = cell_max_[0] - cell_min_[0] + 1;
			int y = cell_max_[1] - cell_min_[1] + 1;
			int z = cell_max_[2] - cell_min_[2] + 1;
			int cell_num = x * y * z;
			dt_ = dt_max_;
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
			/*
			std::cout << "Master Cells\n";
			for (auto it = master_->cellBegin(); it != master_->cellEnd(); it++)
				std::cout << master_->getCellNum(it) << "\n\n";
			std::cout << "pos min: " << domain_cfg_.getCellNum(contact_pmin_) << std::endl;
			std::cout << "pos max: " << domain_cfg_.getCellNum(contact_pmax_) << std::endl;
			std::cout << "end of master Cells\n";
			*/
			for (int i = 0; i < cell_num; i++) {
				Vec3i cur(i % x, (i % (x * y)) / x, i / (x * y));
				Vec3i cell = cur + cell_min_;
				/*
				std::cout << "cell_min_:\n" << cell_min_ << std::endl;
				std::cout << "current cells\n";
				std::cout << cell << std::endl << std::endl;
				*/
				if (master_->hasCell(cell)) {
					auto it_m = master_->findCellStart(cell);
					auto vectors = slave_->getAllAdjacentCells(cell);
					for (int i = master_->begin(it_m); i < master_->end(it_m); i++) {
						if (!vectors.empty()) {
							for (auto sc_it : vectors) {
								for (int j = slave_->begin(sc_it); j < slave_->end(sc_it); j++) {
									(this->*hptr)(i, j);
								}
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
	/// Update contact zone
	//==============================================================================
	void ContactManager::updateContactZone() {
		ensemble_s_->calcBox();
		ensemble_m_->calcBox();
		Vec3d s_pmin = ensemble_s_->getPosMin();
		Vec3d s_pmax = ensemble_s_->getPosMax();
		Vec3d m_pmin = ensemble_m_->getPosMin();
		Vec3d m_pmax = ensemble_m_->getPosMax();

		Vec3i cell_offset = domain_cfg_.getCellOffset();
		Vec3d aug_s_pmin = s_pmin - 1.0 * epsilon_ * cell_offset.cast<double>();
		Vec3d aug_s_pmax = s_pmax + 1.0 * epsilon_ * cell_offset.cast<double>();
		//contact_pmin_ = m_pmin - epsilon_ * cell_offset.cast<double>();
		//contact_pmax_ = m_pmax + epsilon_ * cell_offset.cast<double>();
		contact_pmin_ = aug_s_pmin.cwiseMax(m_pmin);
		contact_pmax_ = aug_s_pmax.cwiseMin(m_pmax);

		cell_min_ = domain_cfg_.getCellNum(contact_pmin_);
		cell_max_ = domain_cfg_.getCellNum(contact_pmax_);
		
	}

	//==============================================================================
	/// Penalty handler
	//==============================================================================
	void ContactManager::penaltyHandler(int i, int j) {
		if ((spos_[j] - mpos_[i]).norm() <  epsilon_ ) {
			Vec3d eta = spos_[j] - mpos_[i];
			Vec3d dvel = svel_[j] - mvel_[i];
			if (dvel.dot(eta) < 0) {
				dt_ = std::min(eta.squaredNorm() / abs(dvel.dot(eta)), dt_);
			}
			Vec3d force = force_const_ * dv_ * eta / eta.norm() * pow( epsilon_ / (eta.norm()) - 1, 1.5);
			macc_[i] -= force / ensemble_m_->getMass();
			sacc_[j] += force / ensemble_s_->getMass();
			//
		}
	}

}