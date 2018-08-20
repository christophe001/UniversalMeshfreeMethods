/*! \file m_model_vtk.h */

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
#include "m_ensemble_creator.h"

namespace msl {
	class ModelVtk : public MsObj {
	protected:
		std::string proj_;
		std::string targ_;
		Vec3d		p_dims_;
		Vec3d		t_dims_;
		double		dp_;
		double		disp_;
		std::string folder_;
		std::string file_;
	public:
		ModelVtk();
		~ModelVtk() {}
		void config(std::string proj, Vec3d p, std::string targ, Vec3d t, double dp, double disp);
		void configIO(std::string folder, std::string file);
		void saveVtk();
	};
}
