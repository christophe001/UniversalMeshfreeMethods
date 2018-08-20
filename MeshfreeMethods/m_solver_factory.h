/*! \file m_solver_factory.h */

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

#include "m_meshfree_base.h"

namespace msl {
	class SolverFactory :public MsObj {
	public:
		~SolverFactory();
		static std::vector<AttributeInfo> getEnsembleAttrInfos(std::string);
		static std::vector<std::string>	getBondAttrNames(std::string);
	};
}
