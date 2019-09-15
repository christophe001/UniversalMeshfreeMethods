/*! \file m_lagrangian_compute.cpp */

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

#include "m_lagrangian_compute.h"
#include <omp.h>

namespace msl {
	LagrangianCompute::LagrangianCompute()
		: ensemble_ptr_(0), sorted_ptr_(0), domain_cfg_(0),
		pos_(0), vel_(0), acc_(0), np_(0),
		cpn_(0), pnbh_(0), nbh_(0) 
	{
		class_name_ = "LagrangianCompute";
	}

	LagrangianCompute::LagrangianCompute(std::shared_ptr<ComputeNeighbor> cpn_ptr)
		: LagrangianCompute::LagrangianCompute()
	{
		init(cpn_ptr);
		class_name_ = "LagrangianCompute";
	}

	LagrangianCompute::LagrangianCompute(std::shared_ptr<SortEnsemble> sorted_ptr)
		: LagrangianCompute::LagrangianCompute()
	{
		init(sorted_ptr);
		class_name_ = "LagrangianCompute";
	}

	void LagrangianCompute::init(std::shared_ptr<SortEnsemble> sorted_ptr) {
		ensemble_ptr_ = sorted_ptr->getEnsemble();
		sorted_ptr_ = sorted_ptr;
		domain_cfg_ = sorted_ptr->domain_cfg_;
		pos_ = ensemble_ptr_->getPos();
		vel_ = ensemble_ptr_->getVel();
		acc_ = ensemble_ptr_->getAcc();
		dict_ = ensemble_ptr_->getDict();
		id_ = ensemble_ptr_->getId();
		np_ = ensemble_ptr_->getNp();
	}

	void LagrangianCompute::init(std::shared_ptr<ComputeNeighbor> cpn_ptr) {
		ensemble_ptr_ = cpn_ptr->ensemble_ptr_;
		sorted_ptr_ = cpn_ptr->sorted_ptr_;
		domain_cfg_ = cpn_ptr->domain_cfg_;
		cpn_ = cpn_ptr;
		pos_ = ensemble_ptr_->getPos();
		vel_ = ensemble_ptr_->getVel();
		acc_ = ensemble_ptr_->getAcc();
		np_ = ensemble_ptr_->getNp();
		pnbh_ = cpn_ptr->computed_ ? cpn_ptr->getPeriNeighborData() : 0;
		nbh_ = cpn_ptr->computed_ ? cpn_ptr->getNeighborhoodData() : 0;
	}

	void LagrangianCompute::computeBond() {
		cpn_ = std::make_shared<ComputeNeighbor>(sorted_ptr_);
		cpn_->compute();
		nbh_ = cpn_->getNeighborhoodData();
		pnbh_ = cpn_->getPeriNeighborData();
	}

	void LagrangianCompute::enforceNoSlip(const Vec3d & center, const Vec3d & vel, 
		const Vec3d & acc, double eps) {
		if (damage_ == nullptr && !ensemble_ptr_->hasScalarAttribute("damage"))
			throwException("enforceNoSlip", "No damage variable found!");
		damage_ = ensemble_ptr_->getScalarAttrPtr("damage")->getAttr();
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
		for (int i = 0; i < np_; i++) {
			Vec3d eta = center - pos_[dict_[i]];
			if (damage_[i] > 0.99 && eta.norm() < eps) {
				acc_[dict_[i]] = acc;
				vel_[dict_[i]] = vel;
			}
		}
	}

	void LagrangianCompute::compute(fpn fun1, fpnpb fun2) {
		if (domain_cfg_ == 0)
			throwException("compute", "fatal error");
		computeAll(fun1);
	}

	void LagrangianCompute::computeAll(fpn funptr) {
		if (nbh_ == 0)
			throwException("computeAll", "No neighborhood data found");
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++) {
			for (auto it = nbh_->begin(i); it != nbh_->end(i); it++) {
				(this->*funptr)(i, it);
			}
		}
	}

	void LagrangianCompute::computeAllwithPeri(fpn fun1, fpnpb fun2) {
		if (nbh_ == 0 || pnbh_ == 0)
			throwException("computeAllWithPeri", "No neighborhood data/perineighbor data found");
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++) {
			for (auto it = nbh_->begin(i); it != nbh_->end(i); it++) {
				(this->*fun1)(i, it);
			}
		}
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < pnbh_->getTotalBond(); i++) {
			PeriNeighborData::PeriBond& pb = pnbh_->getPeriBondList()[i];
			(this->*fun2)(pb);
		}
	}
}