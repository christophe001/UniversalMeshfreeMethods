/*! \file m_self_contact_manager.h */

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

#ifndef _M4_SELF_CONTACT_MANAGER_
#define _M4_SELF_CONTACT_MANAGER_

#include "m_eulerian_compute.h"

namespace msl {
	class SelfContactManager : public EulerianCompute {
	protected:
		Vec3d*	ipos_;
		double	epsilon_;
		double	range_;
		double	ds_max_;
		double	force_const_;
		double	dv_;

	public:
		SelfContactManager(std::shared_ptr<Ensemble> emsemble_ptr, double force_const);
		
		void setParams(Vec3d pmin, Vec3d pmax, double range, 
			double epsilon, double dp);
		
		void setParams(DomainConfig cfg, double epsilon, double dp);
		
		void calcDsMax();
		
		void updateVerletList();
		
		void contactForce(int i, int j);
		
		void computeForces() override;
		
		virtual ~SelfContactManager();
	};
}



#endif // !_M4_SELF_CONTACT_MANAGER_
