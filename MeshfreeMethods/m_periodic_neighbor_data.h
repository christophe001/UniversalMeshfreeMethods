/*! \file m_periodic_neighbor_data.h */

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

#ifndef _M4_PERI_NEIGHBOR_DATA_
#define _M4_PERI_NEIGHBOR_DATA_

#include "m_obj.h"
#include <unordered_map>
#include "m_ensemble_attribute.h"
#pragma warning (disable:4244)

namespace msl {

	class PeriNeighborData : protected MsObj {

	public:
		typedef std::shared_ptr<BondAttribute> BondAttrPtr;

		struct PeriBond {
			int m_i;
			int m_j;
			Vec3d m_offset;
			long m_order;
			PeriBond() : m_i(-1), m_j(-1), m_offset(Vec3d::Zero()), m_order(0) {}
			PeriBond(int i, int j, Vec3d v, long order) 
				: m_i(i), m_j(j), m_offset(v), m_order(order) {}
			//EIGEN_MAKE_ALIGNED_OPERATOR_NEW
			PeriBond& operator=(const PeriBond& other) {
				m_i = other.m_i;
				m_j = other.m_j;
				m_offset = other.m_offset;
				m_order = other.m_order;
				return *this;
			}
		};

		PeriNeighborData()
			: pb_list_(0), list_size_(0), 
			bond_damage_(0), total_bonds_(0)
		{
			class_name_ = "PeriNeighborData";
		}

		PeriNeighborData& operator=(const PeriNeighborData& rhs) = delete;

		~PeriNeighborData() {
			if (pb_list_ != 0) delete[] pb_list_;
			if (bond_damage_ != 0) delete[] bond_damage_;
		}

		void setListSize(long ls);

		void addAttr(std::string s = "");

		bool hasPeriBond(int i);

		bool hasStartPosition(int i);

		long begin(int i);

		long end(int i);

		PeriBond* getPeriBondList() { return pb_list_; }

		int getBondCount(int i) { return bond_count_[i]; }

		double getParicleDamage(int i) { return particle_damage_[i]; }

		double* getBondDamage() { return bond_damage_; }
		
		long getTotalBond() const { return total_bonds_; }

		double getMemoryInMB() const;

		BondAttrPtr getBondAttr(std::string s);

		friend class ComputeNeighbor;

	private:
		std::unordered_map<int, long>			start_positions_;
		std::unordered_map<int, int>			size_;
		std::unordered_map<int, int>			bond_count_;
		std::unordered_map<int, double>			particle_damage_;
		PeriBond*								pb_list_;
		double*									bond_damage_;
		long									list_size_;
		long									total_bonds_;
		std::vector<BondAttrPtr>				bond_attrs_;
		std::unordered_map<std::string, int>	code_;
	};

}

#endif // !_M4_PERI_NEIGHBOR_DATA_