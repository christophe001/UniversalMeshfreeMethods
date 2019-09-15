/*! \file m_lagrangian_compute.h */

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

#ifndef _M4_LAGRANGIAN_
#define _M4_LAGRANGIAN_


#include "m_ensemble.h"
#include "m_sort_ensemble.h"
#include "m_neighborhood_data.h"
#include "m_periodic_neighbor_data.h"
#include "m_compute_neighborhood.h"

namespace msl {
	class LagrangianCompute : public MsObj {
	protected:
		std::shared_ptr<Ensemble>			ensemble_ptr_;
		std::shared_ptr<SortEnsemble>		sorted_ptr_;
		DomainConfig*						domain_cfg_;
		Vec3d*								pos_;
		Vec3d*								vel_;
		Vec3d*								acc_;
		double*								damage_;
		int*								dict_;
		int*								id_;
		int									np_;
		std::shared_ptr<PeriNeighborData>	pnbh_;
		std::shared_ptr<NeighborhoodData>	nbh_;
		std::shared_ptr<ComputeNeighbor>	cpn_;
	
		typedef	void(LagrangianCompute::*fpn)(int i, long j);
		typedef void(LagrangianCompute::*fpnpb)(PeriNeighborData::PeriBond& pb);
		void computeAll(fpn funptr);
		void computeAllwithPeri(fpn fun1, fpnpb fun2);


	public:
		void doNothing(PeriNeighborData::PeriBond& pb) {}
		std::shared_ptr<Ensemble> getEnsemble() { return ensemble_ptr_; }
		LagrangianCompute();
		LagrangianCompute(std::shared_ptr<ComputeNeighbor> cpn_ptr);
		LagrangianCompute(std::shared_ptr<SortEnsemble> sorted_ptr);
		virtual void configModelParams(const std::vector<std::vector<double>>& params) {}
		virtual std::string format() const { return ""; }
		void init(std::shared_ptr<SortEnsemble> sorted_ptr);
		void init(std::shared_ptr<ComputeNeighbor> cpn_ptr);
		void computeBond();
		void enforceNoSlip(const Vec3d& center, const Vec3d& vel, const Vec3d& acc, double eps);
		virtual void compute(fpn fun1, fpnpb fun2 = &doNothing);
		virtual void computeForces() {}
		virtual ~LagrangianCompute() {}
	};
}


#endif // !_M4_LAGRANGIAN_