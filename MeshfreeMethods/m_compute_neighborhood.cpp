/*! \file m_compute_neighborhood.cpp */

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

#include "m_compute_neighborhood.h"
#include <iostream>

namespace msl {

	//==============================================================================
	/// Ctor given ensemble ptr
	//==============================================================================
	ComputeNeighbor::ComputeNeighbor(std::shared_ptr<Ensemble> ensemble_ptr)
		: ensemble_ptr_(ensemble_ptr), np_(ensemble_ptr->getNp()),
		domain_cfg_(0), horizon_(0.0), nb_list_size_(0), pb_count_(0),
		pb_start_(0), pb_list_size_(0), total_cells_(0)
	{
		np_ = ensemble_ptr_->getNp();
		pos_ = ensemble_ptr_->getPos();
		vel_ = ensemble_ptr_->getVel();
		acc_ = ensemble_ptr_->getAcc();
		sorted_ptr_ = std::make_shared<SortEnsemble>(ensemble_ptr);
		nbh_ = std::make_shared<NeighborhoodData>();
		pnbh_ = std::make_shared<PeriNeighborData>();
		computed_ = false;
		class_name_ = "ComputeNeighbor";
	}

	//==============================================================================
	/// Ctor given sorted ptr
	//==============================================================================
	ComputeNeighbor::ComputeNeighbor(std::shared_ptr<SortEnsemble> sorted_ptr)
		: sorted_ptr_(sorted_ptr), ensemble_ptr_(sorted_ptr->getEnsemble()), horizon_(0.0),
		nb_list_size_(0), pb_list_size_(0), pb_count_(0), pb_start_(0)
	{
		np_ = sorted_ptr_->getNp();
		total_cells_ = sorted_ptr_->getTotalCells();
		pos_ = ensemble_ptr_->getPos();
		vel_ = ensemble_ptr_->getVel();
		acc_ = ensemble_ptr_->getAcc();
		domain_cfg_ = new DomainConfig(*sorted_ptr->getDomainConfig());
		nbh_ = std::make_shared<NeighborhoodData>();
		pnbh_ = std::make_shared<PeriNeighborData>();
		computed_ = false;
		class_name_ = "ComputeNeighbor";
		if (domain_cfg_->periActive()) {
			try {
				pb_count_ = new int[np_]();
				pb_start_ = new long[np_]();
			}
			catch (std::bad_alloc) {
				throwException("ComputeNeighbor", "Error occured while allocating memory");
			}
		}

	}

	//==============================================================================
	/// Set domain configuration
	//==============================================================================
	void ComputeNeighbor::setDomainConfig(const DomainConfig& domain_cfg) {
		if (domain_cfg_ != 0)
			delete domain_cfg_;
		domain_cfg_ = new DomainConfig(domain_cfg);
		sorted_ptr_->setDomainConfig(domain_cfg);
		if (domain_cfg_->periActive()) {
			try {
				pb_count_ = new int[np_]();
				pb_start_ = new long[np_]();
			}
			catch (std::bad_alloc) {
				throwException("ComputeNeighbor", "Error occured while allocating memory");
			}
		}
	}

	//==============================================================================
	/// Copy from SortEnsemble
	//==============================================================================
	void ComputeNeighbor::makeSortFull(bool reorder) {
		sorted_ptr_->makeSortFull(reorder);
		total_cells_ = sorted_ptr_->getTotalCells();
	}

	//==============================================================================
	/// Add attribute to bond(scalar)
	//==============================================================================
	void ComputeNeighbor::addAttr(std::string s) {
		nbh_->addAttr(s);
		if (domain_cfg_->periActive())
			pnbh_->addAttr(s);
	}

	//==============================================================================
	/// Initialize meory and stuff
	//==============================================================================
	void ComputeNeighbor::initialize(double dr, double horizon, std::vector<std::string> names) {
		if (!names.empty()) {
			for (auto s : names)
				addAttr(s);
		}
		if (domain_cfg_ == 0)
			throwException("initialize", "no domain configuration found");
		if ((domain_cfg_->getVoxelSize().array() < horizon).any())
			throwException("initialize", "horizon length exceeds cell size");
		horizon_ = horizon;
		bool case2d = domain_cfg_->sim2d();
		int lmax = int(horizon / dr) + 1;
		int ncnt = 0;
		if (case2d) {
			for (int i = -lmax; i <= lmax; i++) {
				for (int j = -lmax; j <= lmax; j++) {
					Vec3d vec{ (double)i, 0.0, (double)j };
					if (vec != Vec3d::Zero() && vec.norm() < horizon / dr) {
						ncnt++;
					}
				}
			}
		}
		else {
			for (int i = -lmax; i <= lmax; i++) {
				for (int j = -lmax; j <= lmax; j++) {
					for (int k = -lmax; k <= lmax; k++) {
						Vec3d vec{ (double)i, (double)k, (double)j };
						if (vec != Vec3d::Zero() && vec.norm() < horizon / dr) {
							ncnt++;
						}
					}
				}
			}
		}
		std::cout << "ncnt: " << ncnt << " np: " << np_ << std::endl;
		long ls = ncnt * np_;
		nb_list_size_ = ls;
		//test
		std::cout << "list size: " << ls << std::endl;
		//test
		nbh_->setNp(np_);
		nbh_->setListSize(nb_list_size_);
	}

	//==============================================================================
	/// Compute neighborhood data and perineighbor data
	//==============================================================================
	void ComputeNeighbor::compute() {
		if (domain_cfg_ == 0)
			throwException("compute", "no domain configuration found");
		if (total_cells_ == 0)
			throwException("compute", "did not sort");
		if (nb_list_size_ == 0)
			throwException("compute", "did not initialize");
		if (domain_cfg_->periActive())
			computeAll();
		else {
			computeNeighbor();
		}
		computed_ = true;
	}

	//==============================================================================
	/// Update particle damage
	//==============================================================================
	void ComputeNeighbor::updateDamage() {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++)
			nbh_->particle_damage_[i] = 0.0;

#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++) {
			for (auto it = nbh_->begin(i); it != nbh_->end(i); it++) {
				nbh_->particle_damage_[i] += nbh_->bond_damage_[it];
			}
		}
		if (pnbh_ != 0) {
			for (auto it : pnbh_->particle_damage_)
				it.second = 0;
			for (auto it : pnbh_->start_positions_) {
				int i = it.first;
				for (long j = pnbh_->begin(i); j != pnbh_->end(i); j++)
					pnbh_->particle_damage_[i] += pnbh_->bond_damage_[j];
			}
		}

	}

	//==============================================================================
	/// Adding split/crack to computed bond structrue
	/// split, cut define split plane cut plane normal direction, sval and cval
	/// defined as val = a1*x + b1*y + c1*z where (a1, b1, c1) is the normal. and
	/// (x, y, z) a point on plane.
	//==============================================================================
	void ComputeNeighbor::addSplit(Vec3d split, double sval, Vec3d cut, double cval) {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++) {
			if (pos_[i].dot(cut) > cval) {
				double s1 = pos_[i].dot(split) - sval;
				for (auto it = nbh_->begin(i); it != nbh_->end(i); it++) {
					int j = nbh_->neighborhood_list_[it];
					double s2 = pos_[j].dot(split) - sval;
					if (s1 * s2 <= 0) {
						nbh_->bond_damage_[it] = 1;
					}
				}
			}
		}
	}

	//==============================================================================
	/// Compute neighborhood
	//==============================================================================
	void ComputeNeighbor::computeNeighbor() {
		int tcells = (sorted_ptr_->getOutNp() == 0) ? total_cells_ : total_cells_ - 1;
		int* nl = nbh_->neighborhood_list_;
		long* sp = nbh_->start_positions_;
		int* bc = nbh_->bond_count_;
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
#endif // _WITH_OMP_
		for (int cid = 0; cid < tcells; cid++) {
			Vec3i cell = sorted_ptr_->getCellNum(cid);
			std::vector<SortEnsemble::CellIt> nb_list = sorted_ptr_->getAllAdjacentCells(cell);
			for (int i = sorted_ptr_->begin(cid); i < sorted_ptr_->end(cid); i++) {
				for (int j = sorted_ptr_->begin(cid); j < sorted_ptr_->end(cid); j++) {
					if (i != j && (pos_[i] - pos_[j]).norm() < 0.9999 * horizon_) {
						bc[i]++;
					}
				}
				if (nb_list.size() != 0) {
					for (auto&& nb : nb_list) {
						for (int j = sorted_ptr_->begin(nb); j < sorted_ptr_->end(nb); j++)
							if ((pos_[i] - pos_[j]).norm() < 0.9999 * horizon_) {
								bc[i]++;
							}
					}
				}
			}
		}
		nbh_->total_bonds_ = 0;
		for (int i = 0; i < np_; i++) {
			sp[i] = nbh_->total_bonds_;
			nbh_->total_bonds_ += bc[i];
		}

#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++)
			bc[i] = 0;

#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
#endif // _WITH_OMP
		for (int cid = 0; cid < tcells; cid++) {
			Vec3i cell = sorted_ptr_->getCellNum(cid);
			std::vector<SortEnsemble::CellIt> nb_list = sorted_ptr_->getAllAdjacentCells(cell);
			for (int i = sorted_ptr_->begin(cid); i < sorted_ptr_->end(cid); i++) {
				for (int j = sorted_ptr_->begin(cid); j < sorted_ptr_->end(cid); j++) {
					if (i != j && (pos_[i] - pos_[j]).norm() < 0.9999 * horizon_) {
						nl[sp[i] + bc[i]] = j;
						bc[i]++;
					}
				}
				if (nb_list.size() != 0) {
					for (auto&& nb : nb_list) {
						for (int j = sorted_ptr_->begin(nb); j < sorted_ptr_->end(nb); j++)
							if ((pos_[i] - pos_[j]).norm() < 0.9999 * horizon_) {
								nl[sp[i] + bc[i]] = j;
								bc[i]++;
							}
					}
				}
			}
		}
	}

	//==============================================================================
	/// Compute neighborhood and perineighborhood
	//==============================================================================
	void ComputeNeighbor::computeAll() {
		computeNeighbor();
		std::cout << "neighborhood data calculated!" << std::endl;
		if (domain_cfg_->periActive()) {
			pb_list_size_ = nb_list_size_ - nbh_->getTotalBonds();
			pnbh_->setListSize(pb_list_size_);
			int tcells = (sorted_ptr_->getOutNp() == 0) ? total_cells_ : total_cells_ - 1;
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
#endif // _WITH_OMP_
			for (int cid = 0; cid < tcells; cid++) {
				Vec3i cell = sorted_ptr_->getCellNum(cid);
				if (domain_cfg_->cellOnSurface(cell)) {
					std::vector<SortEnsemble::PeriCellIt> pb_list = sorted_ptr_->getAllPeriAdjCells(cell);
					if (pb_list.size() != 0) {
						for (int i = sorted_ptr_->begin(cid); i < sorted_ptr_->end(cid); i++) {
							for (auto&& pb : pb_list) {
								SortEnsemble::CellIt pit = pb.first;
								Vec3d offset = pb.second;
								for (int j = sorted_ptr_->begin(pit); j < sorted_ptr_->end(pit); j++) {
									if ((pos_[i] - pos_[j] - offset).norm() < 0.999 * horizon_) {
										pb_count_[i]++;
									}
								}
							}
						}
					}
				}
			}
			pnbh_->total_bonds_ = 0;
			for (int i = 0; i < np_; i++) {
				pb_start_[i] = pnbh_->total_bonds_;
				if (pb_count_[i] != 0) {
					(pnbh_->bond_count_)[i] = pb_count_[i];
					(pnbh_->size_)[i] = pb_count_[i];
					pb_count_[i] = 0;
					(pnbh_->start_positions_)[i] = pnbh_->total_bonds_;
					pnbh_->total_bonds_ += pnbh_->bond_count_[i];
				}
			}
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
#endif // _WITH_OMP_
			for (int cid = 0; cid < tcells; cid++) {
				Vec3i cell = sorted_ptr_->getCellNum(cid);
				if (domain_cfg_->cellOnSurface(cell)) {
					std::vector<SortEnsemble::PeriCellIt> pb_list = sorted_ptr_->getAllPeriAdjCells(cell);
					if (pb_list.size() != 0) {
						for (int i = sorted_ptr_->begin(cid); i < sorted_ptr_->end(cid); i++) {
							for (auto&& pb : pb_list) {
								SortEnsemble::CellIt pit = pb.first;
								Vec3d offset = pb.second;
								for (int j = sorted_ptr_->begin(pit); j < sorted_ptr_->end(pit); j++) {
									if ((pos_[i] - pos_[j] - offset).norm() < 0.999 * horizon_) {
										long order = pb_start_[i] + pb_count_[i];
										pb_count_[i]++;
										PeriNeighborData::PeriBond pb(i, j, offset, order);
										(pnbh_->pb_list_)[order] = pb;
									}
								}
							}
						}
					}
				}
			}
		}
	}
}