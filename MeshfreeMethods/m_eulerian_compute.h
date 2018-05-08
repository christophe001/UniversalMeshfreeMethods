/*! \file m_eulerian_compute.h */

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

#ifndef _M4_EULERIAN_
#define _M4_EULERIAN_

#include "m_ensemble.h"
#include "m_sort_ensemble.h"
#include <functional>

namespace msl {
	class EulerianCompute : public MsObj {
	protected:
		std::shared_ptr<Ensemble>		ensemble_ptr_;
		std::shared_ptr<SortEnsemble>	sorted_ptr_;
		Vec3d*							pos_;
		Vec3d*							vel_;
		Vec3d*							acc_;
		int								np_;
		int*							dict_;
		int*							id_;
		int								total_cells_;

		// interaction functions, specified in derived class
		typedef void(EulerianCompute::*fpn)(int i, int j);
		//typedef std::function<void(int, int)> fpn;

		// compute cell interaction with itself
		void computeSelf(SortEnsemble::CellIt it, fpn funptr);
		void computeSelf(int cid, fpn funptr);

		// compute cell interaction with adjacent cells
		// (with cell_id greater than current cell)
		void computeAdjacent(SortEnsemble::CellIt it, fpn funptr);
		void computeAdjacent(int cid, fpn funptr);


		// compute cell interaction with periodically adjacent cells(if applicable)
		void computePeriAdj(SortEnsemble::CellIt it, fpn funptr);
		void computePeriAdj(int cid, fpn funptr);


		// compute cell interaction with itself and with adjacent cells 
		// (with cell_id greater than current cell)
		void computeAll(SortEnsemble::CellIt it, fpn funptr);
		void computeAll(int cid, fpn funptr);


		// compute all interactions with current cell,
		// including periodically neighbored cells
		void computeAllWithPeri(SortEnsemble::CellIt it, fpn funptr1, fpn funptr2);
		void computeAllWithPeri(int cid, fpn funptr1, fpn funptr2);

	public:
		EulerianCompute(std::shared_ptr<Ensemble> emsemble_ptr);
		void setDomainConfig(DomainConfig domain_cfg);
		virtual ~EulerianCompute() {}

		virtual std::string format() const = 0;

		void makeSortFull(bool reorder);
		void makeSortPartial(bool all = true);
		void doNothing(int i, int j) {}
		virtual void compute(fpn fun1, fpn fun2 = &doNothing);
		virtual void computeForces() {}
	};
}

#endif // !_M4_LAGRANGIAN_