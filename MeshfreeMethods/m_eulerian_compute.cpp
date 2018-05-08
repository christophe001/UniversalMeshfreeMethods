/*! \file m_eulerian_compute.cpp */

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

#include "m_eulerian_compute.h"
#include <omp.h>

namespace msl {
	EulerianCompute::EulerianCompute(std::shared_ptr<Ensemble> ensemble_ptr)
		: ensemble_ptr_(ensemble_ptr), np_(ensemble_ptr->getNp()), total_cells_(0)
	{
		sorted_ptr_ = std::make_shared<SortEnsemble>(ensemble_ptr);
		pos_ = ensemble_ptr_->getPos();
		vel_ = ensemble_ptr_->getVel();
		acc_ = ensemble_ptr_->getAcc();
		dict_ = ensemble_ptr_->getDict();
		id_ = ensemble_ptr_->getId();
		class_name_ = "EulerianCompute";
	}

	void EulerianCompute::setDomainConfig(DomainConfig domain_cfg) {
		sorted_ptr_->setDomainConfig(domain_cfg);
	}

	void EulerianCompute::makeSortFull(bool reorder) {
		sorted_ptr_->makeSortFull(reorder);
		total_cells_ = sorted_ptr_->getTotalCells();
	}

	void EulerianCompute::makeSortPartial(bool all) {
		sorted_ptr_->makeSortPartial(all);
		total_cells_ = sorted_ptr_->getTotalCells();
	}

	void EulerianCompute::computeSelf(SortEnsemble::CellIt it, fpn funptr) {
		for (int i = sorted_ptr_->begin(it); i < sorted_ptr_->end(it) - 1; i++) {
			for (int j = i + 1; j < sorted_ptr_->end(it); j++)
				if (i != j) {
					(this->*funptr)(i, j);
				}
		}
	}

	void EulerianCompute::computeSelf(int cid, fpn funptr) {
		for (int i = sorted_ptr_->begin(cid); i < sorted_ptr_->end(cid) - 1; i++)
			for (int j = i + 1; j < sorted_ptr_->end(cid); j++)
				if (i != j) {
					(this->*funptr)(i, j);
				}
	}

	void EulerianCompute::computeAdjacent(SortEnsemble::CellIt it, fpn funptr) {
		Vec3i cell = sorted_ptr_->getCellNum(it);
		std::vector<SortEnsemble::CellIt> nb_list = sorted_ptr_->getAllAdjacentCells(cell);
		if (nb_list.size() != 0) {
			for (int i = sorted_ptr_->begin(it); i < sorted_ptr_->end(it); i++) {
				for (auto&& nb : nb_list) {
					for (int j = sorted_ptr_->begin(nb); j < sorted_ptr_->end(nb); j++) {
						(this->*funptr)(i, j);
					}
				}
			}
		}
	}

	void EulerianCompute::computeAdjacent(int cid, fpn funptr) {
		Vec3i cell = sorted_ptr_->getCellNum(cid);
		std::vector<SortEnsemble::CellIt> nb_list = sorted_ptr_->getAllAdjacentCells(cell);
		if (nb_list.size() != 0) {
			for (int i = sorted_ptr_->begin(cid); i < sorted_ptr_->end(cid); i++) {
				for (auto&& nb : nb_list) {
					for (int j = sorted_ptr_->begin(nb); j < sorted_ptr_->end(nb); j++) {
						(this->*funptr)(i, j);
					}
				}
			}
		}
	}

	void EulerianCompute::computePeriAdj(SortEnsemble::CellIt it, fpn funptr) {
		Vec3i cell = sorted_ptr_->getCellNum(it);
		std::vector<SortEnsemble::PeriCellIt> pb_list = sorted_ptr_->getAllPeriAdjCells(cell);
		if (pb_list.size() != 0) {
			for (int i = sorted_ptr_->begin(it); i < sorted_ptr_->end(it); i++) {
				for (auto&& pb : pb_list) {
					SortEnsemble::CellIt pit = pb.first;
					// Vec3d offset = pb.second;
					for (int j = sorted_ptr_->begin(pit); j < sorted_ptr_->end(pit); j++)
						(this->*funptr)(i, j);
				}
			}
		}
	}

	void EulerianCompute::computePeriAdj(int cid, fpn funptr) {
		Vec3i cell = sorted_ptr_->getCellNum(cid);
		std::vector<SortEnsemble::PeriCellIt> pb_list = sorted_ptr_->getAllPeriAdjCells(cell);
		if (pb_list.size() != 0) {
			for (int i = sorted_ptr_->begin(cid); i < sorted_ptr_->end(cid); i++) {
				for (auto&& pb : pb_list) {
					SortEnsemble::CellIt pit = pb.first;
					// Vec3d offset = pb.second;
					for (int j = sorted_ptr_->begin(pit); j < sorted_ptr_->end(pit); j++)
						(this->*funptr)(i, j);
				}
			}
		}
	}

	void EulerianCompute::computeAll(SortEnsemble::CellIt it, fpn funptr) {
		Vec3i cell = sorted_ptr_->getCellNum(it);
		std::vector<SortEnsemble::CellIt> nb_list = sorted_ptr_->getAllAdjacentCells(cell);
		for (int i = sorted_ptr_->begin(it); i < sorted_ptr_->end(it); i++) {
			for (int j = i + 1; j < sorted_ptr_->end(it); j++) {
				(this->*funptr)(i, j);
			}
			if (nb_list.size() != 0) {
				for (auto&& nb : nb_list) {
					for (int j = sorted_ptr_->begin(nb); j < sorted_ptr_->end(nb); j++)
						(this->*funptr)(i, j);
				}
			}
		}
	}

	void EulerianCompute::computeAll(int cid, fpn funptr) {
		Vec3i cell = sorted_ptr_->getCellNum(cid);
		std::vector<SortEnsemble::CellIt> nb_list = sorted_ptr_->getAllAdjacentCells(cell);
		for (int i = sorted_ptr_->begin(cid); i < sorted_ptr_->end(cid); i++) {
			for (int j = i + 1; j < sorted_ptr_->end(cid); j++) {
				(this->*funptr)(i, j);
			}
			if (nb_list.size() != 0) {
				for (auto&& nb : nb_list) {
					for (int j = sorted_ptr_->begin(nb); j < sorted_ptr_->end(nb); j++)
						(this->*funptr)(i, j);
				}
			}
		}
	}


	void EulerianCompute::computeAllWithPeri(SortEnsemble::CellIt it, fpn funptr1, fpn funptr2) {
		Vec3i cell = sorted_ptr_->getCellNum(it);
		std::vector<SortEnsemble::CellIt> nb_list = sorted_ptr_->getAllAdjacentCells(cell);
		std::vector<SortEnsemble::PeriCellIt> pb_list = sorted_ptr_->getAllPeriAdjCells(cell);
		for (int i = sorted_ptr_->begin(it); i < sorted_ptr_->end(it); i++) {
			for (int j = i + 1; j < sorted_ptr_->end(it); j++) {
				(this->*funptr1)(i, j);
			}
			if (nb_list.size() != 0) {
				for (auto&& nb : nb_list) {
					for (int j = sorted_ptr_->begin(nb); j < sorted_ptr_->end(nb); j++)
						(this->*funptr1)(i, j);
				}
			}
			if (pb_list.size() != 0) {
				for (auto&& pb : pb_list) {
					SortEnsemble::CellIt pit = pb.first;
					Vec3d offset = pb.second;
					for (int j = sorted_ptr_->begin(pit); j < sorted_ptr_->end(pit); j++)
						(this->*funptr2)(i, j);
				}
			}
		}
	}

	void EulerianCompute::computeAllWithPeri(int cid, fpn funptr1, fpn funptr2) {
		Vec3i cell = sorted_ptr_->getCellNum(cid);
		std::vector<SortEnsemble::CellIt> nb_list = sorted_ptr_->getAllAdjacentCells(cell);
		std::vector<SortEnsemble::PeriCellIt> pb_list = sorted_ptr_->getAllPeriAdjCells(cell);
		for (int i = sorted_ptr_->begin(cid); i < sorted_ptr_->end(cid); i++) {
			for (int j = i + 1; j < sorted_ptr_->end(cid); j++) {
				(this->*funptr1)(i, j);
			}
			if (nb_list.size() != 0) {
				for (auto&& nb : nb_list) {
					for (int j = sorted_ptr_->begin(nb); j < sorted_ptr_->end(nb); j++)
						(this->*funptr1)(i, j);
				}
			}
			if (pb_list.size() != 0) {
				for (auto&& pb : pb_list) {
					SortEnsemble::CellIt pit = pb.first;
					Vec3d offset = pb.second;
					for (int j = sorted_ptr_->begin(pit); j < sorted_ptr_->end(pit); j++)
						(this->*funptr2)(i, j);
				}
			}
		}
	}

	
	void EulerianCompute::compute(fpn fun1, fpn fun2) {
		if (total_cells_ == 0)
			throwException("compute", "did not sort");
		int tcells = (sorted_ptr_->getOutNp() == 0) ? total_cells_ : total_cells_ - 1;
		if (sorted_ptr_->getDomainConfig()->periActive()) {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif
			for (int cid = 0; cid < tcells; cid++) {
				computeAllWithPeri(cid, fun1, fun2);
			}
		}
		else {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif
			for (int cid = 0; cid < tcells; cid++) {
				computeAll(cid, fun1);
			}
		}
	}


}