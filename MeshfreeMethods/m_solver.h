/*! \file m_solver.h */

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

#ifndef _M4_SOLVER_
#define _M4_SOLVER_

#include <iostream>
#include "m_eulerian_compute.h"
#include "m_lagrangian_compute.h"
#include "m_compute_neighborhood.h"

#pragma warning (disable:4244)

namespace msl {
	class ComputeSolver : protected MsObj {
	protected:
		std::shared_ptr<Ensemble>				ensemble_;
		std::shared_ptr<SortEnsemble>			sorted_;
		std::shared_ptr<NeighborhoodData>		nbh_;
		std::shared_ptr<PeriNeighborData>		pbh_;
		std::shared_ptr<ComputeNeighbor>		compute_;
		Vec3d*									pos_;
		Vec3d*									vel_;
		Vec3d*									acc_;
		std::vector<Ensemble::ScalarAttrPtr>	scalar_attrs_;
		std::vector<Ensemble::VectorAttrPtr>	vector_attrs_;
		int										np_;
		double									dp_;
		double									horizon_, dt_, T_;
		std::string								filename_;
		std::string								folder_;
		int										timestep_, max_step_;
		int										sv_step_, sv_count_;

	public:
		ComputeSolver();
		~ComputeSolver() {}
		void printLabel();
		void printInfo();
		void printMemoryInfo();
		void createCase(std::string folder, std::string xml) {}
		void create2DCasePeridynamics(Peridm2DParams pp, bool ec = true);
		double getMaxAcc() const;
		double getMaxVel() const;
		Vec3d getMaxPos(Vec3i dir) const;
		void initialize() {}
		double getModulus();
		void verletUpdate();
		void configParams(std::string folder, std::string filename, double dt, double total_time);
		void saveVtk();
		void run();

	};
}


#endif // !_M4_SOLVER_

