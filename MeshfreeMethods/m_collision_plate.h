/*! \file m_collision_plate.h */

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

#include "m_meta_solver.h"
#include "m_state_kl_model.h"

namespace msl {

	class CollisionPlate : public MsObj {
	protected:
		ObjInfo projectile_;
		ObjInfo target_;
		std::shared_ptr<MetaSolver> meta_;

	public:
		CollisionPlate(ObjInfo projectile,  ObjInfo target, DomainConfig domain_cfg, 
			std::string folder, std::string file, double dt, double total_time, int sv = 30);
		virtual ~CollisionPlate() {}
		void config();
		void run();

	};
}
