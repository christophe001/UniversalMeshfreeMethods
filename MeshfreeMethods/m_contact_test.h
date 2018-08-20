/*! \file m_contact_test.h */

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
#ifndef _M4_CONTACT_TEST_
#define _M4_CONTACT_TEST_

#include "m_collision_plate.h"

namespace msl {
	class ContactTest : public MsObj {
	protected:
		//*********************************************************************
		//		projectile							target
		//			\							____  /
		//				____   displacement    |    |
		//		       |	\	____________   |    |
		//		       |____/			       |    |
		//									   |____|
		//
		//*********************************************************************
		std::string		proj_shape_;
		Vec3d			proj_dims_;
		double			proj_vel_;
		std::string		proj_solv_;

		double			disp_;
		std::string		targ_shape_;
		Vec3d			targ_dims_;
		std::string		targ_solv_;

		double			dp_;
		double			dt_, total_time_;
		int				sv_steps_;
		
		std::string		folder_;
		std::string		file_;

		std::shared_ptr<MetaSolver> solver_;
		int				config_state_;
	
	public:
		ContactTest();
		~ContactTest() {}
		void configShapes(std::string proj, Vec3d p_dims, std::string targ, Vec3d t_dims);
		void configSolver(std::string proj = "rigid", std::string targ = "kl");
		void configParams(double vel, double disp, double dp,std::string folder, std::string file);
		void configTime(double dt, double total_time, int sv);
		void initialize();
		void run();

	};
}


#endif
