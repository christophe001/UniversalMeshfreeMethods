/*! \file m_sort_ensemble.cpp */

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

#include "m_sort_ensemble.h"
#include <omp.h>
#include <algorithm>
#include "u_timer.h"

using msl::mortonEncode;

namespace msl {
	//==============================================================================
	/// Ctor
	//==============================================================================
	SortEnsemble::SortEnsemble(std::shared_ptr<Ensemble> ensemble_ptr)
		: ensemble_ptr_(ensemble_ptr), np_(ensemble_ptr->getNp()),
		start_positions_(0), total_cells_(0),
		out_np_(0), valid_np_(0), domain_cfg_(0)
	{
		class_name_ = "SortEnsemble";
		pos_ = ensemble_ptr_->getPos();
		vel_ = ensemble_ptr_->getVel();
		acc_ = ensemble_ptr_->getAcc();
		id_ = ensemble_ptr_->getId();
		dict_ = ensemble_ptr_->getDict();
		scalar_attrs_ = ensemble_ptr_->getScalarAttrPtrs();
		vector_attrs_ = ensemble_ptr_->getVectorAttrPtrs();
		tensor_attrs_ = ensemble_ptr_->getTensorAttrPtrs();
		if (np_) {
			try {
				cells_ = new uint64_t[np_];
				new_order_ = new int[np_];
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
				for (int i = 0; i < np_; i++)
					new_order_[i] = i;

				int sz = tensor_attrs_.size() == 0 ? sizeof(Vec3d) : sizeof(Mat3d);
				long mem = long(sz) * long(np_);
				sort_stack_ = new unsigned char[mem];
			}
			catch (std::bad_alloc) {
				throwException("Constructor", "Error occured while allocating memory");
			}
		}

	}

	//==============================================================================
	/// Destructor
	//==============================================================================
	SortEnsemble::~SortEnsemble() {
		if (domain_cfg_ != 0 )		delete domain_cfg_;
		if (start_positions_ != 0)	delete[] start_positions_;
		if (cells_ != 0)			delete[] cells_;
		if (new_order_ != 0)		delete[] new_order_;
		if (sort_stack_ != 0)		delete[] sort_stack_;
	}

	//==============================================================================
	/// Set domain configuration from another DomainConfig object
	//==============================================================================
	void SortEnsemble::setDomainConfig(DomainConfig cfg) {
		if (domain_cfg_ != 0)
			delete domain_cfg_;
		domain_cfg_ = new DomainConfig(cfg);
	}

	//==============================================================================
	/// Translate out domain periodic particles so that all particles in domain
	//==============================================================================
	void SortEnsemble::makePeriInDomain() {
		if (domain_cfg_ == 0) return;
		if (domain_cfg_->periActive()) {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif
			for (int i = 0; i < np_; i++) {
				domain_cfg_->periTranslate(pos_[i]);
			}
		}
	}

	//==============================================================================
	/// Reorder particle id, dict is used to indicate position of particle with id
	/// i in current sorted ensemble
	//==============================================================================
	void SortEnsemble::reorderID() {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++) {
			id_[i] = i;
			dict_[i] = i;		//! given particle id, find its offset in memory
		}
	}

	//==============================================================================
	/// Sort cells
	//==============================================================================
	void SortEnsemble::sortCells() {
		// return if no DomainConfig object is set
		if (domain_cfg_ == 0) return;

		// total cells occupied and out_np_ reset to 0
		total_cells_ = 0;
		out_np_ = 0;

		// deallocate memory for sorting cells
		if (start_positions_ != 0)
			delete[] start_positions_;

		if (!parts_in_cell_.empty())
			parts_in_cell_.clear();

		if (!cell_start_.empty())
			cell_start_.clear();
		
#ifdef _WITH_OMP_
		const int nt = omp_get_max_threads();
		auto pic = new std::unordered_map<uint64_t, int>[nt];
		int out = 0;
#pragma omp parallel for schedule(static) num_threads(nt) reduction(+:out)
		for (int i = 0; i < np_; i++) {
			new_order_[i] = i;
			if (domain_cfg_->isOutBound(pos_[i])) out++;
			cells_[i] = mortonEncode(domain_cfg_->getCellNum(pos_[i]));
			pic[omp_get_thread_num()][cells_[i]]++;
		}
		out_np_ = out;
		valid_np_ = np_ - out_np_;

		for (int i = 0; i < nt; i++) {
			for (auto it = pic[i].begin(); it != pic[i].end(); it++)
				parts_in_cell_[it->first] += it->second;
		}

		total_cells_ = parts_in_cell_.size();
		delete[] pic;
#else
		for (int i = 0; i < np_; i++) {
			new_order_[i] = i;
			if (domain_cfg_->isOutBound(pos_[i]))
				out_np_++;
			cells_[i] = mortonEncode(domain_cfg_->getCellNum(pos_[i]));
			parts_in_cell_[cells_[i]]++;
		}
		valid_np_ = np_ - out_np_;
		total_cells_ = parts_in_cell_.size();
#endif

		try {
			start_positions_ = new int[total_cells_];
		}
		catch (std::bad_alloc) {
			throwException("sortCells", "Error occurred while allocating memory");
		}
		int count = 0, idx = 0;
		for (auto it = parts_in_cell_.begin(); it != parts_in_cell_.end(); it++) {
			cell_start_[it->first] = count;
			start_positions_[idx++] = count;
			count += it->second;
		}
	}

	//==============================================================================
	/// Template method to sort all data types
	//==============================================================================
	template <class T>
	void SortEnsemble::sort(T* vec) {
#ifdef _WITH_OMP_
#pragma omp parallel
{
#pragma omp for schedule(static) 
			for (int i = 0; i < np_; i++)
				((T*)sort_stack_)[i] = vec[new_order_[i]];
#pragma omp for schedule(static)
			for (int j = 0; j < np_; j++)
				vec[j] = ((T*)sort_stack_)[j];
}
#else 
		for (int i = 0; i < np_; i++)
			((T*)sort_stack_)[i] = vec[new_order_[i]];
		for (int j = 0; j < np_; j++)
			vec[j] = ((T*)sort_stack_)[j];
#endif // _WITH_OMP_
	}

	//==============================================================================
	/// Sort ensemble attributes
	//==============================================================================
	void SortEnsemble::sortAttributes() {
		for (auto&& sa : scalar_attrs_)
			sort<double>(sa->getAttr());

		for (auto&& va : vector_attrs_)
			sort<Vec3d>(va->getAttr());

		for (auto&& ta : tensor_attrs_)
			sort<Mat3d>(ta->getAttr());
	}

	//==============================================================================
	/// Make sort full
	//==============================================================================
	void SortEnsemble::makeSortFull(bool reorder) {
		ensemble_ptr_->calcBox();
		sortCells();
		sortParticles();
		sortAttributes();
		if (reorder) reorderID();
	}

	//==============================================================================
	/// Make sort partial
	//==============================================================================
	void SortEnsemble::makeSortPartial(bool all) {
		ensemble_ptr_->calcBox();
		sortCells();
		sortParticles(all);
	}

	//==============================================================================
	/// Sort particles given sorted cells
	//==============================================================================
	void SortEnsemble::sortParticles(bool all) {
		if (cell_start_.empty()) return;
		std::sort(new_order_, new_order_ + np_,
			[this](int& i, int& j) {return cells_[i] < cells_[j]; });
		sort<Vec3d>(pos_);
		if (all) {
			sort<Vec3d>(acc_);
			sort<Vec3d>(vel_);
		}
		sort<int>(id_);
		//sortAttributes();
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif //  _WITH_OMP_
		for (int i = 0; i < np_; i++)
			dict_[id_[i]] = i;
	}

	//==============================================================================
	/// Get adjacent cellits that has larger encoded cell number
	//==============================================================================
	std::vector<SortEnsemble::CellIt> SortEnsemble::getAdjacentCells(Vec3i cell) {
		std::vector<SortEnsemble::CellIt> res;
		if (cell == domain_cfg_->getMapSize() + vec3one)
			return res;
		uint64_t code = mortonEncode(cell);
		Vec3i temp;
		if (domain_cfg_->sim2d()) {
			for (int x = -1; x < 2; x++) {
				for (int z = -1; z < 2; z++) {
					Vec3i inc( x, 0, z );
					if (inc != Vec3i::Zero()) {
						if ((inc.array() <= arr3zero).all()) continue;
						if ((inc.array() >= arr3zero).all()) {
							temp = cell + inc;
							CellIt it = cell_start_.find(mortonEncode(temp));
							if (it != cell_start_.end())
								res.push_back(it);
						}
						else {
							temp = cell + inc;
							uint64_t code_r = mortonEncode(temp);
							if ((temp.array() >= arr3zero).all() && code_r > code) {
								CellIt it = cell_start_.find(code_r);
								if (it != cell_start_.end())
									res.push_back(it);
							}
						}
					}
				}
			}
		}
		else {
			for (int x = -1; x < 2; x++) {
				for (int y = -1; y < 2; y++) {
					for (int z = -1; z < 2; z++) {
						Vec3i inc( x, y, z );
						if (inc != Vec3i::Zero()) {
							if ((inc.array() <= arr3zero).all()) continue;
							if ((inc.array() >= arr3zero).all()) {
								temp = cell + inc;
								CellIt it = cell_start_.find(mortonEncode(temp));
								if (it != cell_start_.end())
									res.push_back(it);
							}
							else {
								temp = cell + inc;
								uint64_t code_r = mortonEncode(temp);
								if ((temp.array() >= arr3zero).all() && code_r > code) {
									CellIt it = cell_start_.find(code_r);
									if (it != cell_start_.end())
										res.push_back(it);
								}
							}
						}
					}
				}
			}
		}
		return res;
	}

	//==============================================================================
	/// Get all adjacent cells given current cell
	//==============================================================================
	std::vector<SortEnsemble::CellIt> SortEnsemble::getAllAdjacentCells(Vec3i cell) {
		std::vector<SortEnsemble::CellIt> res;
		uint64_t code = mortonEncode(cell);
		if (domain_cfg_->sim2d()) {
			for (int x = -1; x < 2; x++) {
				for (int z = -1; z < 2; z++) {
					Vec3i inc( x, 0, z );
					if (inc != Vec3i::Zero()) {
						Vec3i temp = cell + inc;
						if ((temp.array() >= arr3zero).all()) {
							CellIt it = cell_start_.find(mortonEncode(temp));
							if (it != cell_start_.end())
								res.push_back(it);
						}
					}
				}
			}
		}
		else {
			for (int x = -1; x < 2; x++) {
				for (int y = -1; y < 2; y++) {
					for (int z = -1; z < 2; z++) {
						Vec3i inc( x, y, z );
						if (inc != Vec3i::Zero()) {
							Vec3i temp = cell + inc;
							if ((temp.array() >= arr3zero).all()) {
								CellIt it = cell_start_.find(mortonEncode(temp));
								if (it != cell_start_.end())
									res.push_back(it);
							}
						}
					}
				}
			}
		}
		return res;
	}

	//==============================================================================
	/// Get periodically adjacent cells that has larger encoded cell number
	//==============================================================================
	std::vector<SortEnsemble::PeriCellIt> SortEnsemble::getPeriAdjCells(Vec3i cell) {
		std::vector<SortEnsemble::PeriCellIt> pclist;
		Vec3d domain_size = domain_cfg_->getDomainSize();
		if (domain_cfg_->cellOnSurface(cell)) {
			if (domain_cfg_->sim2d()) {
				for (int i = -1; i <= 1; i++) {
					for (int j = -1; j <= 1; j++) {
						Vec3i offset = Vec3i(i, 0, j);
						if (offset != Vec3i::Zero()) {
							Vec3i temp = cell + offset;
							if (domain_cfg_->cellInPeri(temp)) {
								Vec3i peri = temp.cwiseMin(Vec3i::Zero()) +
									(temp - domain_cfg_->getMaxCell()).cwiseMax(Vec3i::Zero());
								Vec3i pcell = temp - peri.cwiseProduct(domain_cfg_->getMapSize());
								uint64_t code = mortonEncode(pcell);
								if (code > mortonEncode(cell)) {
									CellIt it = cell_start_.find(mortonEncode(pcell));
									if (it != cell_start_.end())
										pclist.push_back(std::make_pair(it,
											domain_size.cwiseProduct(peri.cast<double>())));
								}
							}
						}
					}
				}
			}
			else {
				for (int i = -1; i <= 1; i++) {
					for (int j = -1; j <= 1; j++) {
						for (int k = -1; k <= 1; k++) {
							Vec3i offset = Vec3i(i, j, k);
							if (offset != Vec3i::Zero()) {
								Vec3i temp = cell + offset;
								if (domain_cfg_->cellInPeri(temp)) {
									Vec3i peri = temp.cwiseMin(Vec3i::Zero()) +
										(temp - domain_cfg_->getMaxCell()).cwiseMax(Vec3i::Zero());
									Vec3i pcell = temp - peri.cwiseProduct(domain_cfg_->getMapSize());
									uint64_t code = mortonEncode(pcell);
									if (code > mortonEncode(cell)) {
										CellIt it = cell_start_.find(mortonEncode(pcell));
										if (it != cell_start_.end())
											pclist.push_back(std::make_pair(it,
												domain_size.cwiseProduct(peri.cast<double>())));
									}
								}
							}
						}
					}
				}
			}
		}
		return pclist;
	}
	

	//==============================================================================
	/// Get all periodically adjacent cells
	//==============================================================================
	std::vector<SortEnsemble::PeriCellIt> SortEnsemble::getAllPeriAdjCells(Vec3i cell) {
		std::vector<SortEnsemble::PeriCellIt> pclist;
		Vec3d domain_size = domain_cfg_->getDomainSize();
		if (domain_cfg_->cellOnSurface(cell)) {
			if (domain_cfg_->sim2d()) {
				for (int i = -1; i <= 1; i++) {
					for (int j = -1; j <= 1; j++) {
						Vec3i offset = Vec3i(i, 0, j);
						if (offset != Vec3i::Zero()) {
							Vec3i temp = cell + offset;
							if (domain_cfg_->cellInPeri(temp)) {
								Vec3i peri = temp.cwiseMin(Vec3i::Zero()) +
									(temp - domain_cfg_->getMaxCell()).cwiseMax(Vec3i::Zero());
								Vec3i pcell = temp - peri.cwiseProduct(domain_cfg_->getMapSize());
								CellIt it = cell_start_.find(mortonEncode(pcell));
								if (it != cell_start_.end())
									pclist.push_back(std::make_pair(it, 
									domain_size.cwiseProduct(peri.cast<double>())));
							}
						}
					}
				}
			}
			else {
				for (int i = -1; i <= 1; i++) {
					for (int j = -1; j <= 1; j++) {
						for (int k = -1; k <= 1; k++) {
							Vec3i offset = Vec3i(i, j, k);
							if (offset != Vec3i::Zero()) {
								Vec3i temp = cell + offset;
								if (domain_cfg_->cellInPeri(temp)) {
									Vec3i peri = temp.cwiseMin(Vec3i::Zero()) +
										(temp - domain_cfg_->getMaxCell()).cwiseMax(Vec3i::Zero());
									Vec3i pcell = temp - peri.cwiseProduct(domain_cfg_->getMapSize());
									CellIt it = cell_start_.find(mortonEncode(pcell));
									if (it != cell_start_.end())
										pclist.push_back(std::make_pair(it, 
											domain_size.cwiseProduct(peri.cast<double>())));
								}
							}
						}
					}
				}
			}
		}
		return pclist;
	}


	//==============================================================================
	/// Get memory used for SortEnsemble in MB
	//==============================================================================
	double SortEnsemble::getMemoryInMB() const {
		int64_t size = sizeof(int) * (2 * np_ + 3) + np_ * sizeof(Vec3d) * sizeof(unsigned char)
			+ np_ * sizeof(uint64_t) + 4 * sizeof(int*) + 3 * sizeof(Vec3d*) + sizeof(DomainConfig*)
			+ sizeof(uint64_t*) + sizeof(unsigned char*);
		double memsize = size / 1048576.0;
		return memsize;
	}

}
