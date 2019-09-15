/*! \file m_compute_neighborhood.h */

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

#ifndef _M4_COMPUTE_NEIGHBOR_
#define _M4_COMPUTE_NEIGHBOR_

#include "m_ensemble.h"
#include "m_sort_ensemble.h"
#include "m_neighborhood_data.h"
#include "m_periodic_neighbor_data.h"

namespace msl {
	class ComputeNeighbor : protected MsObj {
	protected:
		std::shared_ptr<Ensemble>			ensemble_ptr_;
		std::shared_ptr<SortEnsemble>		sorted_ptr_;
		DomainConfig*						domain_cfg_;
		Vec3d*								pos_;
		Vec3d*								vel_;
		Vec3d*								acc_;
		int									np_;
		long								total_cells_;
		std::shared_ptr<PeriNeighborData>	pnbh_;
		std::shared_ptr<NeighborhoodData>	nbh_;
		double								horizon_;
		long								nb_list_size_;
		long								pb_list_size_;
		bool								computed_;
		int*								pb_count_;
		long*								pb_start_;

		// compute neighborhood data
		void computeNeighbor();

		// compute neighborhood data & periodic neighborhood data 
		void computeAll();

	public:
		friend class LagrangianCompute;
		friend class EnsembleCreator;

		ComputeNeighbor(std::shared_ptr<Ensemble> ensemble_ptr);
		
		ComputeNeighbor(const ComputeNeighbor& cn) = delete;
		
		ComputeNeighbor& operator=(const ComputeNeighbor& cn) = delete;

		ComputeNeighbor(std::shared_ptr<SortEnsemble> sorted_ptr);
		
		~ComputeNeighbor() {
			if (domain_cfg_ != 0)
				delete domain_cfg_;
			if (pb_count_ != 0)
				delete[] pb_count_;
			if (pb_start_ != 0)
				delete[] pb_start_;
		}
		void setDomainConfig(const DomainConfig& domain_cfg);
		
		void addAttr(std::string s = "");

		void makeSortFull(bool reorder = true);

		void initialize(double dr, double horizon, std::vector<std::string> names);

		void compute();

		void updateDamage();

		void addSplit(Vec3d split, double sval, Vec3d cut, double cval);	//! adding split/crack to computed bond structure.

		std::shared_ptr<NeighborhoodData> getNeighborhoodData() { return nbh_; }

		std::shared_ptr<PeriNeighborData> getPeriNeighborData() { return pnbh_; }

	};
}

#endif