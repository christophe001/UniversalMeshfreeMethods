/*! \file m_contact_test.cpp */

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
#include "m_contact_test.h"

namespace msl {
	ContactTest::ContactTest() {
		config_state_ = 0;
		class_name_ = "ContactTest";
	}

	void ContactTest::configShapes(std::string proj, Vec3d p_dims, std::string targ, Vec3d t_dims) {
		proj_shape_ = proj;
		proj_dims_ = p_dims;
		targ_shape_ = targ;
		targ_dims_ = t_dims;
		config_state_ |= 1;
	}

	void ContactTest::configSolver(std::string proj, std::string targ) {
		proj_solv_ = proj;
		targ_solv_ = targ;
		config_state_ |= 2;
	}

	void ContactTest::configParams(double vel, double disp, double dp, std::string folder, std::string file) {
		proj_vel_ = vel;
		disp_ = disp;
		dp_ = dp;
		folder_ = folder;
		file_ = file;
		config_state_ |= 4;
	}

	void ContactTest::configTime(double dt, double total_time, int sv) {
		dt_ = dt;
		total_time_ = total_time;
		sv_steps_ = sv;
		config_state_ |= 8;
	}

	void ContactTest::initialize() {
		if (config_state_ != 15)
			throwException("initialize", "not all parameters configured");
		if (solver_ != nullptr)
			solver_ = nullptr;
		ObjInfo proj("sphere", proj_dims_, "rigid", dp_, 8000, Vec3d{ 0, 0, -proj_vel_ }, Vec3d{0, 0, disp_});
		ObjInfo targ("rectangle", targ_dims_, "kl", dp_, 2200);
		solver_ = std::make_shared<CollisionPlate>(proj, targ);
	}

	void ContactTest::run() {
		initialize();
		solver_->run();
	}
}
