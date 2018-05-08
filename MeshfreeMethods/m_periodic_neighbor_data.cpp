/*! \file m_periodic_neighbor_data.cpp */

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

#include "m_periodic_neighbor_data.h"

namespace msl {
	void PeriNeighborData::setListSize(long ls)
	{
		if (pb_list_ != 0)
			delete[] pb_list_;
		if (bond_damage_ != 0)
			delete[] bond_damage_;
		try {
			pb_list_ = new PeriBond[ls];
			bond_damage_ = new double[ls];
			list_size_ = ls;
		}
		catch (std::bad_alloc) {
			throwException("setListSize()", "Error occured while allocating memory");
		}
		for (auto && attr : bond_attrs_)
			attr->setListSize(ls);
	}

	void PeriNeighborData::addAttr(std::string s) {
		int sz = code_.size();
		auto attr_ptr = std::make_shared<BondAttribute>();
		bond_attrs_.push_back(std::move(attr_ptr));
		code_[s] = sz;
	}

	bool PeriNeighborData::hasPeriBond(int i) {
		return (bond_count_.find(i) != bond_count_.end()) && (bond_count_[i] != 0);
	}

	bool PeriNeighborData::hasStartPosition(int i) {
		return (start_positions_.find(i) != start_positions_.end());
	}

	long PeriNeighborData::begin(int i) {
#ifdef _DEBUG_
		if (!hasPeriBond(i))
			throwException("end()", "particle doesn't have peri bond");
#endif
		return start_positions_[i];
	}

	long PeriNeighborData::end(int i) {
#ifdef _DEBUG_
		if (!hasPeriBond(i))
			throwException("end()", "particle doesn't have peri bond");
#endif
		return start_positions_[i] + bond_count_[i];
	}

	double PeriNeighborData::getMemoryInMB() const {
		uint64_t size = list_size_ * sizeof(PeriBond) + sizeof(PeriBond*) + 3 * sizeof(int);
		double memsize = size / 1048576.0;
		for (auto&& ba : bond_attrs_) {
			memsize += ba->getMemoryInMB();
		}
		return memsize;
	}

	PeriNeighborData::BondAttrPtr PeriNeighborData::getBondAttr(std::string s) {
		if (code_.find(s) == code_.end())
			throwException("getBondAttr", "Cannot find attribute");
		return bond_attrs_[code_[s]];
	}
}