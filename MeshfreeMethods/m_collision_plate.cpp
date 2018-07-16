/*! \file m_collision_plate.cpp */

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

#include "m_collision_plate.h"

namespace msl {
	CollisionPlate::CollisionPlate(ObjInfo projectile, ObjInfo target, 
		std::vector<AttributeInfo>& p_infos, std::vector<AttributeInfo>& t_infos, DomainConfig domain_cfg,
		double dt, double total_time, int sv)
		: projectile_(projectile), target_(target) {
		setDomainConfig(domain_cfg);
		configRun(dt, total_time, sv);
		class_name_ = "CollisionPlate";
		createScene(p_infos, t_infos);
	}

	void CollisionPlate::createScene(std::vector<AttributeInfo>& p_infos,
		std::vector<AttributeInfo>& t_infos) {
		if (projectile_.dp != target_.dp)
			throwException("createScene", "Unmatched dp!");
		dp_ = projectile_.dp;
		slave_solv_->setDomainConfig(domain_config_);
		master_solv_->setDomainConfig(domain_config_);
		slave_solv_->createEnsemble(projectile_, domain_config_, p_infos);
		master_solv_->createEnsemble(target_, domain_config_, t_infos);
		/*
		s_creator_ = std::make_shared<EnsembleCreator>(projectile_.shape, projectile_.offset);
		s_creator_->setDims(projectile_.dims);
		s_creator_->setOrientation(projectile_.orientation, projectile_.theta);
		s_creator_->setDpDensity(projectile_.dp, projectile_.density);

		m_creator_ = std::make_shared<EnsembleCreator>(target_.shape, target_.offset);
		m_creator_->setDims(target_.dims);
		m_creator_->setOrientation(target_.orientation, target_.theta);
		m_creator_->setDpDensity(target_.dp, target_.density);
		*/
		configMeshfreeSolver();
		configContact();
	}

	void CollisionPlate::configMeshfreeSolver() {
		std::vector<std::string> emp;
		MetaSolver::configMeshfreeSolver("rigid", "", emp, "kl", "", emp);
	}

	void CollisionPlate::configContact() {
		MetaSolver::configContact();
		cm_ = std::make_shared<ContactManager>(master_sort_, slave_sort_, dp_, dt_max_);
	}
}