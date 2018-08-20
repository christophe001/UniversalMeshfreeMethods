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
	CollisionPlate::CollisionPlate(ObjInfo projectile, ObjInfo target,  DomainConfig domain_cfg,
		std::string folder, std::string file, double dt, double total_time, int sv)
		: projectile_(projectile), target_(target) {
		meta_ = std::make_shared<MetaSolver>();
		meta_->configIO(folder, file);
		meta_->configRun(dt, total_time, sv);
		meta_->createScene(projectile_, target_, domain_cfg);
		meta_->configContact();
		class_name_ = "CollisionPlate";
	}

	void CollisionPlate::run() {
		meta_->run();
	}

}