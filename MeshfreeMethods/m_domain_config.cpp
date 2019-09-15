/*! \file m_domain_config.cpp */

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

#include "m_domain_config.h"

namespace msl {

	//==============================================================================
	/// Ctor
	//==============================================================================
	DomainConfig::DomainConfig(const Vec3d& pos_min, const Vec3d& pos_max,
		const double& cell_size, const PeriType& peri)
		: pmin_(pos_min), pmax_(pos_max), domain_size_(pos_max - pos_min)
	{
		class_name_ = "DomainConfig";
		if ((pmin_.array() > pmax_.array()).any())
			throwException("Constructor", "Error in setting limits of domain");
		if (cell_size <= 0)
			throwException("Constructor", "Invalid cell size, must > 0");

		case2d_ = (pmax_[1] == pmin_[1]);
		int p = static_cast<int>(peri);
		peribc_ = Vec3i( p & 1, (p & 2) >> 1, (p & 4) >> 2 );
		Vec3d cell = 1.005 * cell_size * Vec3d::Ones();
		cell_offset_ = vec3one;
		map_size_ = domain_size_.cwiseQuotient(cell).cast<int>();
		if (case2d_) { map_size_[1] = 1; cell_offset_[1] = 0; }
		voxel_size_ = domain_size_.cwiseQuotient(map_size_.cast<double>());
		if (case2d_) voxel_size_[1] = cell_size;
		peri_offset_ = peribc_.cwiseProduct(cell_offset_);
		peri_min_ = pmin_ - voxel_size_.cwiseProduct(peri_offset_.cast<double>());
		peri_max_ = pmax_ + voxel_size_.cwiseProduct(peri_offset_.cast<double>());
		max_cell_ = map_size_ - vec3one;
		peri_cellmin_ = vec3zero - peri_offset_;
		peri_cellmax_ = max_cell_ + peri_offset_;
	}

	//==============================================================================
	/// Copy constructor
	//==============================================================================
	DomainConfig::DomainConfig(const DomainConfig& other)
		: pmin_(other.pmin_),
		pmax_(other.pmax_),
		peri_min_(other.peri_min_),
		peri_max_(other.peri_max_),
		domain_size_(other.domain_size_),
		voxel_size_(other.voxel_size_),
		map_size_(other.map_size_),
		case2d_(other.case2d_),
		cell_offset_(other.cell_offset_),
		peri_offset_(other.peri_offset_),
		peribc_(other.peribc_),
		max_cell_(other.max_cell_),
		peri_cellmin_(other.peri_cellmin_),
		peri_cellmax_(other.peri_cellmax_)
	{
		class_name_ = "DomainConfig";
	}

	//==============================================================================
	/// Copy operator = 
	//==============================================================================
	DomainConfig& DomainConfig::operator= (const DomainConfig& other) {
		this->pmin_ = other.pmin_;
		this->pmax_ = other.pmax_;
		this->peri_min_ = other.peri_min_;
		this->peri_max_ = other.peri_max_;
		this->domain_size_ = other.domain_size_;
		this->voxel_size_ = other.voxel_size_;
		this->map_size_ = other.map_size_;
		this->case2d_ = other.case2d_;
		this->cell_offset_ = other.cell_offset_;
		this->peri_offset_ = other.peri_offset_;
		this->peribc_ = other.peribc_;
		this->max_cell_ = other.max_cell_;
		this->peri_cellmin_ = other.peri_cellmin_;
		this->peri_cellmax_ = other.peri_cellmax_;
		this->class_name_ = other.getClassName();
		return *this;
	}

	//==============================================================================
	/// bool operator ==
	//==============================================================================
	bool DomainConfig::operator==(const DomainConfig & other) {
		if (this->pmin_ != other.pmin_) return false;
		if (this->pmax_ != other.pmax_) return false;
		if (this->peri_min_ != other.peri_min_) return false;
		if (this->peri_max_ != other.peri_max_) return false;
		if (this->domain_size_ != other.domain_size_) return false;
		if (this->voxel_size_ != other.voxel_size_) return false;
		if (this->map_size_ != other.map_size_) return false;
		if (this->case2d_ != other.case2d_) return false;
		if (this->cell_offset_ != other.cell_offset_) return false;
		if (this->peri_offset_ != other.peri_offset_) return false;
		if (this->peribc_ != other.peribc_) return false;
		if (this->max_cell_ != other.max_cell_) return false;
		if (this->peri_cellmin_ != other.peri_cellmin_) return false;
		if (this->peri_cellmax_ != other.peri_cellmax_) return false;
		return true;
	}

	//==============================================================================
	/// Set cell size
	//==============================================================================
	void DomainConfig::setCellSize(double s) {
		Vec3d cell = 1.005 * s * Vec3d::Ones();
		map_size_ = domain_size_.cwiseQuotient(cell).cast<int>();
		cell_offset_ = vec3one;
		if (case2d_) { map_size_[1] = 1; cell_offset_[1] = 0; }
		voxel_size_ = domain_size_.cwiseQuotient(map_size_.cast<double>());
		if (case2d_) voxel_size_[1] = s;
		peri_offset_ = peribc_.cwiseProduct(cell_offset_);
		peri_min_ = pmin_ - voxel_size_.cwiseProduct(peri_offset_.cast<double>());
		peri_max_ = pmax_ + voxel_size_.cwiseProduct(peri_offset_.cast<double>());
		max_cell_ = map_size_ - vec3one;
		peri_cellmin_ = vec3zero - peri_offset_;
		peri_cellmax_ = max_cell_ + peri_offset_;
	}

	//==============================================================================
	/// Get cell number given position
	//==============================================================================
	Vec3i DomainConfig::getCellNum(const Vec3d& pos) const {
		if (isOutBound(pos))
			return map_size_ + vec3one;
		Vec3d rpos = pos - peri_min_;
		return rpos.cwiseQuotient(voxel_size_).cast<int>() - peri_offset_;
	}

	//==============================================================================
	/// Function peri translate
	//==============================================================================
	void DomainConfig::periTranslate(Vec3d& pos) const {
		if (!inPeriZone(pos)) return;
		Vec3i cell = getCellNum(pos);
		Vec3i cmin = vec3zero.cwiseMin(cell);
		Vec3i tp = map_size_ - vec3one;
		Vec3i cmax = tp.cwiseMax(cell) - map_size_ + vec3one;
		pos -= voxel_size_.cwiseProduct(map_size_.cast<double>()).cwiseProduct(cmin.cast<double>() + cmax.cast<double>());
	}
	
	//==============================================================================
	/// ostream << operator, display info
	//==============================================================================
	std::ostream& operator<<(std::ostream& os, const DomainConfig& cfg) {
		os << "Domain range: \n"
			<< "Min: " << cfg.getPosMin().transpose().format(CleanFmt) 
			<< "Max: " << cfg.getPosMax().transpose().format(CleanFmt) << "\n"
			<< "Periodic range:\n"
			<< "Min: " << cfg.getPeriMin().transpose().format(CleanFmt) 
			<< "Max: " << cfg.getPeriMax().transpose().format(CleanFmt) << "\n"
			<< "( Conditions:  1, True   0, False )\n"
			<< "Sim2D: " << cfg.sim2d() << "\n"
			<< "PeriX: " << cfg.getPeriOffset()[0]
			<< "\t PeriY: " << cfg.getPeriOffset()[1]
			<< "\t PeriZ: " << cfg.getPeriOffset()[2]
			<< std::endl;

		os << "Size of domain: " << cfg.getDomainSize().transpose().format(CleanFmt)
			<< "\nVoxel size: " << cfg.getVoxelSize().transpose().format(CleanFmt)
			<< "\nCell map size: " << cfg.getMapSize().transpose().format(CleanFmt) << std::endl;
		return os;
	}
}