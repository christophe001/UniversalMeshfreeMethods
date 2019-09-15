/*! \file m_sort_ensemble.h */

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

#ifndef _M4_SORT_ENSEMBLE_
#define _M4_SORT_ENSEMBLE_

#include "m_ensemble.h"
#include "m_functions.h"
#include "m_domain_config.h"
#include <memory>
#include <map>
#include <unordered_map>
#pragma warning (disable:4244)

namespace msl {
	class SortEnsemble : protected MsObj {
	protected:
		int												np_;
		int												valid_np_;
		int												out_np_;
		long											total_cells_;
		std::shared_ptr<Ensemble>						ensemble_ptr_;

		DomainConfig*									domain_cfg_;
		std::map<uint64_t, int>							parts_in_cell_;
		std::unordered_map<uint64_t, int>				cell_start_;
		Vec3d*											pos_;
		Vec3d*											vel_;
		Vec3d*											acc_;
		std::vector<std::shared_ptr<ScalarAttribute> >	scalar_attrs_;
		std::vector<std::shared_ptr<VectorAttribute> >  vector_attrs_;
		std::vector<std::shared_ptr<TensorAttribute> >  tensor_attrs_;
		int*											id_;				//! particle id
		int*											dict_;				//! current particle # given particle id

		int*											start_positions_;
		int*											new_order_;
		unsigned char*									sort_stack_;
		uint64_t*										cells_;

	public:
		friend class LagrangianCompute;
		
		friend class ContactManager;
		
		SortEnsemble(std::shared_ptr<Ensemble> ensemble_ptr);
		
		virtual ~SortEnsemble();

		void setDomainConfig(DomainConfig cfg);
		
		void makePeriInDomain();
		
		void reorderID();

		virtual void sortCells();
		
		virtual void sortParticles(bool all = true);

		int getNp() const { return np_; }
		
		int getTotalCells() const { return total_cells_; }
		
		int getOutNp() const { return out_np_; }
		
		int getValidNp() const { return valid_np_; }


		DomainConfig* getDomainConfig() { return domain_cfg_; }
		
		std::shared_ptr<Ensemble> getEnsemble() { return ensemble_ptr_; }

		template<class T> void sort(T* vec);
		void sortAttributes();

		void makeSortFull(bool reorder);
		
		void makeSortPartial(bool all = false);		//! Sort only pos vel acc, and doesn't reorder ids
		
		typedef std::unordered_map<uint64_t, int>::iterator CellIt;
		
		typedef std::pair<CellIt, Vec3d> PeriCellIt;
		
		CellIt cellBegin() { return cell_start_.begin(); }
		
		CellIt cellEnd() { return cell_start_.end(); }
		
		CellIt findCellStart(Vec3i cellnum) { return cell_start_.find(mortonEncode(cellnum)); }
		
		bool hasCell(Vec3i cellnum) { return cell_start_.find(mortonEncode(cellnum)) != cell_start_.end(); }
		
		CellIt findCellStart(uint64_t code) { return cell_start_.find(code); }
		
		bool hasCell(uint64_t code) { return cell_start_.find(code) != cell_start_.end(); }
		
		Vec3i getCellNum(CellIt it) { return domain_cfg_->getCellNum(pos_[begin(it)]); }
		
		Vec3i getCellNum(int cell_id) { return domain_cfg_->getCellNum(pos_[start_positions_[cell_id]]); }
		
		std::vector<CellIt> getAdjacentCells(Vec3i cell);
		
		std::vector<CellIt> getAllAdjacentCells(Vec3i cell);
		
		std::vector<PeriCellIt> getPeriAdjCells(Vec3i cell);
		
		std::vector<PeriCellIt> getAllPeriAdjCells(Vec3i cell);
		
		double getMemoryInMB() const;
		
		int begin(CellIt it) { return it->second; }
		
		int end(CellIt it) { return it->second + parts_in_cell_[it->first]; }
		
		int begin(long cell_id) { return start_positions_[cell_id]; }
		
		int end(long cell_id) { return cell_id == total_cells_ - 1 ? np_ : start_positions_[cell_id + 1]; }
	};
}

#endif // !_M4_SORT_ENSEMBLE