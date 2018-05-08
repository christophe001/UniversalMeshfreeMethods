/*! \file m_domain_config.h */

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

#ifndef _M4_DOMAIN_CONFIG_
#define _M4_DOMAIN_CONFIG_

#include "m_types.h"
#include "m_obj.h"
#pragma warning (disable:4244)
#pragma warning (disable:4267)

namespace msl {
	class DomainConfig : protected MsObj {
	private:
		Vec3d     pmin_, pmax_;					//! Define domain range
		Vec3d     peri_min_, peri_max_;         //! Domain range recalculated to offset periodic condition 
		Vec3d     domain_size_;                 //! Size of the domain
		Vec3d     voxel_size_;					//! Size of each voxel cell
		Vec3i     map_size_;                    //! Total voxel map for the domain
		bool      case2d_;                      //! If case is 2d simulation
		Vec3i     cell_offset_, peri_offset_;   //! Cell offsets 
		Vec3i     peribc_;                      //! Periodic bcs
		Vec3i     max_cell_;					//! max cell for particles in domain
		Vec3i     peri_cellmin_, peri_cellmax_; //! defines cells bound for peri zone

	public:
		DomainConfig() : peribc_(vec3zero) { class_name_ = "DomainConfig"; }
		DomainConfig(const Vec3d& pos_min, const Vec3d& pos_max, 
			const double& cell_size, const PeriType& peri = PeriType::kPeriO);
		DomainConfig(const DomainConfig& other);
		DomainConfig& operator= (const DomainConfig& other);
		~DomainConfig() {}
		bool operator==(const DomainConfig& other);
		bool operator!=(const DomainConfig& other) { return !(*this == other); }

		void setCellSize(double s);
		//Vec3u GetCellNum(const Vec3d& pos) const;
		Vec3d getPosMin() const { return pmin_; }
		Vec3d getPosMax() const { return pmax_; }
		Vec3d getPeriMin() const { return peri_min_; }
		Vec3d getPeriMax() const { return peri_max_; }
		Vec3d getDomainSize() const { return domain_size_; }
		Vec3d getVoxelSize() const { return voxel_size_; }
		Vec3i getMapSize() const { return map_size_; }
		Vec3i getMaxCell() const { return max_cell_; }
		bool  sim2d() const { return case2d_; }
		Vec3i getPeriType() const { return peribc_; }
		Vec3i getCellOffset() const { return cell_offset_; }
		Vec3i getPeriOffset() const { return peri_offset_; }
		bool periActive() const { return peribc_ != vec3zero; }
		void periTranslate(Vec3d& pos) const;
		Vec3i getCellNum(const Vec3d& pos) const;
		bool cellOutBound(const Vec3i& cell) const {
			return !((cell.array() >= arr3zero).all() && (cell.array() <= max_cell_.array()).all());
		}
		bool cellInPeri(const Vec3i& cell) const {
			return cellOutBound(cell) && (cell.array() >= peri_cellmin_.array()).all() && 
				(cell.array() <= peri_cellmax_.array()).all();
		}
		bool cellInnerDomain(const Vec3i& cell) const {
			return (cell.array() >= (vec3zero + peri_offset_).array()).all() && 
				(cell.array() <= (max_cell_ - peri_offset_).array()).all();
		}
		bool cellOnSurface(const Vec3i& cell) const {
			return !cellInnerDomain(cell) && !cellOutBound(cell);
		}
		bool isOutBound(const Vec3d& pos) const { return !((pos.array() >= peri_min_.array()).all() 
			&& (pos.array() <= peri_max_.array()).all()); }
		bool isOutDomain(const Vec3d& pos) const { return !((pos.array() >= pmin_.array()).all() 
			&& (pos.array() <= pmax_.array()).all()); }
		bool inPeriZone(const Vec3d& pos) const { return isOutDomain(pos) && (!isOutBound(pos)); }
		friend std::ostream& operator<< (std::ostream& os, const DomainConfig& cfg);
	};
}

#endif

