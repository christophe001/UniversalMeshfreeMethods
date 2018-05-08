/*! \file m_neighborhood_data.h */

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

#ifndef _M4_NEIGHBORHOOD_DATA_
#define _M4_NEIGHBORHOOD_DATA_

#include <memory>
#include <unordered_map>
#include <vector>
#include "m_obj.h"
#include "omp.h"
#include "m_ensemble_attribute.h"
#pragma warning (disable:4244)

namespace msl {

	class NeighborhoodData : protected MsObj {
	public:
		typedef std::shared_ptr<BondAttribute> BondAttrPtr;

		NeighborhoodData()
			: np_(0), list_size_(0), total_bonds_(0), neighborhood_list_(0),
			start_positions_(0), bond_count_(0), bond_damage_(0), particle_damage_(0)
		{
			class_name_ = "NeighborhoodData";
		}

		NeighborhoodData& operator=(const NeighborhoodData& rhs) = delete;

		~NeighborhoodData();

		void setNp(int np);
		
		void setListSize(long nl);

		long begin(int i);

		long end(int i);

		void addAttr(std::string s = "");

		double getMemoryInMB() const;

		int  getNp() const { return np_; }
		
		int  getListSize() const { return list_size_; }
		
		long getTotalBonds() const { return total_bonds_; }
		
		int* getNeighborhoodList() { return neighborhood_list_; }

		double* getBondDamage() { return bond_damage_; }

		double* getParticleDamage() { return particle_damage_; }
		
		int* getBondCount() { return bond_count_; }
		
		BondAttrPtr getBondAttr(std::string s);

		friend class ComputeNeighbor;

	protected:
		int		np_;								//! number of particles
		long	list_size_;							//! list size of neighborhood list
		long	total_bonds_;						//! total bonds
		int*	neighborhood_list_;					//! neighborhood list store information of bonds
		long*	start_positions_;					//! return position in nbh of particle's first bond
		int*	bond_count_;						//! bond count of each particle
		double*	bond_damage_;						//! store damage of each bond, in same order as start positions
		double*	particle_damage_;					//! store total sum of bond damage
		std::vector<BondAttrPtr> bond_attrs_;		//! store bond attribute ptrs
		std::unordered_map<std::string, int> code_; //! map bond attribute name to its location.
	};
}
#endif // !_M4_NEIGHBORHOOD_DATA