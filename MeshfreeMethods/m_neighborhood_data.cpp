/*! \file m_neighborhood_data.cpp */

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

#include "m_neighborhood_data.h"

namespace msl {
	void NeighborhoodData::setListSize(long ls) {
		if (neighborhood_list_ != 0)
			delete[] neighborhood_list_;
		if (bond_damage_ != 0)
			delete[] bond_damage_;
		try {
			neighborhood_list_ = new int[ls];
			bond_damage_ = new double[ls]();
			list_size_ = ls;
		}
		catch (const std::bad_alloc&) {
			throwException("setListSize()", "Error occured while allocating memory");
		}
		for (auto && attr : bond_attrs_)
			attr->setListSize(ls);
	}

	long NeighborhoodData::begin(int i) {
#ifdef _DEBUG_
		if (id < 0 || id >= np_)
			throwException("begin()", "Invalid paritcle id");
#endif
		return start_positions_[i];
	}

	long NeighborhoodData::end(int i) {
#ifdef _DEBUG_
		if (id < 0 || id >= np_)
			throwException("begin()", "Invalid paritcle id");
#endif
		return start_positions_[i] + bond_count_[i];
	}

	void NeighborhoodData::addAttr(std::string s) {
		int sz = bond_attrs_.size();
		auto attr_ptr = std::make_shared<BondAttribute>(s);
		bond_attrs_.push_back(std::move(attr_ptr));
		code_[s] = sz;
	}

	double NeighborhoodData::getMemoryInMB() const {
		int64_t size = sizeof(int) * (3 + list_size_ + np_) + sizeof(long) * np_ 
			+ sizeof(double) * list_size_ + sizeof(double) * np_;
		double memsize = size / 1048576.0;
		for (auto && attr : bond_attrs_)
			memsize += attr->getMemoryInMB();
		return memsize;
	}

	NeighborhoodData::BondAttrPtr NeighborhoodData::getBondAttr(std::string s)
	{
		if (code_.find(s) == code_.end())
			throwException("getBondAttr", "Cannot find attribute");
		return bond_attrs_[code_[s]];
	}

	void NeighborhoodData::setNp(int np) {
		if (start_positions_ != 0)
			delete[] start_positions_;
		if (bond_count_ != 0)
			delete[] bond_count_;
		if (particle_damage_ != 0)
			delete[] particle_damage_;
		try {
			np_ = np;
			start_positions_ = new long[np + 1];
			bond_count_ = new int[np]();
			particle_damage_ = new double[np]();
		}
		catch (const std::bad_alloc&) {
			throwException("setNp()", "Error occured while allocating memory");
		}
	}

	NeighborhoodData::~NeighborhoodData() {
		if (neighborhood_list_ != 0)
			delete[] neighborhood_list_;
		if (start_positions_ != 0)
			delete[] start_positions_;
		if (bond_count_ != 0)
			delete[] bond_count_;
		if (bond_damage_ != 0)
			delete[] bond_damage_;
		if (particle_damage_ != 0)
			delete[] particle_damage_;
	}

}