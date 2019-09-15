/*! \file m_self_contact_manager.cpp */

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

#include "m_self_contact_manager.h"
#include <algorithm>
#include <cmath>

namespace msl {
	//==============================================================================
	/// Ctor
	//==============================================================================
	SelfContactManager::SelfContactManager(std::shared_ptr<Ensemble> emsemble_ptr, double force_const)
		: EulerianCompute(ensemble_ptr_), force_const_(force_const) {
		class_name_ = "SelfContactManager";
		try {
			last_pos_ = new Vec3d[np_];
		}
		catch (std::bad_alloc) {
			throwException("Constructor", "Error occured while allocating memory");
		}
		ds_max_ = 0.0;
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++)
			last_pos_[i] = pos_[i];
	}

	//==============================================================================
	/// Set params for contact
	/// Object creation: ctor -> setParams
	/// Usage:	updateVerletList -> computeRepulsiveForce
	//==============================================================================
	void SelfContactManager::setParams(Vec3d pmin, Vec3d pmax, 
		double range, double epsilon, double dp) {
		DomainConfig cfg(pmin, pmax, range);
		EulerianCompute::setDomainConfig(cfg);
		range_ = range;
		epsilon_ = epsilon;
		dv_ = cfg.sim2d() ? dp * dp : dp * dp * dp;
		makeSortFull(false);
	}

	//==============================================================================
	/// Set params for contact
	/// Object creation: ctor -> setParams
	/// Usage:	updateVerletList -> computeRepulsiveForce
	//==============================================================================
	void SelfContactManager::setParams(DomainConfig cfg, double epsilon, double dp) {
		EulerianCompute::setDomainConfig(cfg);
		range_ = cfg.getVoxelSize()[0];
		epsilon_ = epsilon;
		dv_ = cfg.sim2d() ? dp * dp : dp * dp * dp;
		makeSortFull(false);
	}

	//==============================================================================
	/// Destructor
	//==============================================================================
	SelfContactManager::~SelfContactManager() {
		if (last_pos_ != 0)
			delete[] last_pos_;
	}

	//==============================================================================
	/// Calculate maximum displacement since last updating of verlet list
	//==============================================================================
	void SelfContactManager::calcDsMax() {
		ds_max_ = 0.0;
		for (int i = 0; i < np_; i++) {
			ds_max_ = std::max((pos_[i]- last_pos_[i]).norm(), ds_max_);
		}
	}

	//==============================================================================
	/// If maximum displacement greater than threshold, update particle order
	/// otherwise update maximum displacement
	//==============================================================================
	void SelfContactManager::updateVerletList() {
		calcDsMax();
		if (range_ - epsilon_ - 2.0 * ds_max_ < 0) {
			sorted_ptr_->makeSortPartial(true);
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
			for (int i = 0; i < np_; i++)
				last_pos_[i] = pos_[i];
		}
	}

	//==============================================================================
	/// Pairwise contact force
	//==============================================================================
	void SelfContactManager::contactForce(int i, int j) {
		if ((pos_[j] - pos_[i]).norm() < epsilon_) {
			Vec3d eta = pos_[j] - pos_[i];
			//acc_[i] -= force_const_ * dv_ * eta / (eta.norm() * ensemble_ptr_->getMass()) 
			//	* pow(epsilon_ / (eta.norm()) - 1, 1.0);
			acc_[i] -= force_const_ * dv_ * eta / (eta.norm() * ensemble_ptr_->getMass()) * log(epsilon_ / eta.norm());
		}
	}

	//==============================================================================
	/// compute repulsive(penaulty) force for all particles
	//==============================================================================
	void SelfContactManager::computeForces() {
		updateVerletList();
		compute(static_cast<EulerianCompute::fpn>(&SelfContactManager::contactForce));
	}

	//==============================================================================
	/// formatting the manager
	//==============================================================================
	std::string SelfContactManager::format() const {
		std::string s;
		s += " range is: " + std::to_string(range_) + "\n";
		s += " epsilon is: " + std::to_string(epsilon_) + "\n";
		s += " force constant is: " + std::to_string(force_const_) + "\n";
		return s;
	}
}