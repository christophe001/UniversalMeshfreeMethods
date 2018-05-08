/*! \file m_neighborhood_test.h */

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

#ifndef _M4_NEIGHBORHOOD_TEST_
#define _M4_NEIGHBORHOOD_TEST_

#include "m_compute_neighborhood.h"
#include "m_ensemble_creator.h"

namespace msl {
	class NeighborhoodTest : public MsObj {
	protected:
		bool								sim2d_;
		Vec3d								dims_;
		double								dp_, horizon_;
		std::shared_ptr<EnsembleCreator>	creator_;
		std::shared_ptr<Ensemble>			ensemble_;
		//std::shared_ptr<SortEnsemble>		sorted_;
		std::shared_ptr<ComputeNeighbor>	cpn_;
		std::shared_ptr<NeighborhoodData>	nbh_;
		std::shared_ptr<PeriNeighborData>	pbh_;
		std::string							folder_, file_;
		int									np_;
		Vec3d*								pos_; 
		Vec3d*								vel_;
		Vec3d*								acc_;

	public:
		NeighborhoodTest(bool sim2d, Vec3d dims, double dp);
		virtual ~NeighborhoodTest() {}
		void setHorizon(double horizon) { horizon_ = horizon; }
		void run();	//! create ensemble and compute neighbor
		void runPeri();
		void setSaveParams(std::string folder, std::string file);
		void saveVtk();
	};
}

#endif // !_M4_NEIGHBORHOOD_TEST_

