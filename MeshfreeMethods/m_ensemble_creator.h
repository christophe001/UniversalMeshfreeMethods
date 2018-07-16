/*! \file m_ensemble_creator.h */

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

#ifndef _M4_ENSEMBLE_CREATOR_
#define _M4_ENSEMBLE_CREATOR_
#include "m_ensemble.h"
#include "m_shape.h"

namespace msl {
	class EnsembleCreator : public MsObj {
		std::shared_ptr<Ensemble>	ensemble_;
		std::shared_ptr<Shape>		shape_;
		std::shared_ptr<Shape>		shape_actual_;
		double						dp_;
		double						density_;
		Vec3d						offset_;
		Vec3d*						cache_;
	public:
		EnsembleCreator(std::string name, Vec3d offset = Vec3d::Zero());
		std::shared_ptr<Ensemble>	getEnsemble() {
			return ensemble_;
		}
		EnsembleCreator(const EnsembleCreator& ec) = delete;
		EnsembleCreator& operator=(const EnsembleCreator& ec) = delete;
		virtual ~EnsembleCreator();
		void addScalarAttributes(std::vector<std::string> scalars);
		void addVectorAttributes(std::vector<std::string> vectors);
		void addTensorAttributes(std::vector<std::string> tensors);
		std::shared_ptr<Shape> getActualShape() { return shape_actual_; }
		void setDims(double a, double b = 0, double c = 0) { shape_->setDims(a, b, c); }
		void setDims(Vec3d dims) { 
			if (shape_->shape2D()) shape_->setDims(dims[0], dims[2]); 
			else shape_->setDims(dims[0], dims[1], dims[2]);
		}
		void setOrientation(const Vec3d& dir = Vec3d::UnitZ(), const double& theta = 0) { 
			shape_->setOrientation(dir, theta); 
		}
		void setScalarAttribute(std::string s, double* sa);
		void setScalarAttributeConstant(std::string s, double a);
		void setVectorAttribute(std::string s, Vec3d*  va);
		void setVectorAttributeZero(std::string s);
		void setTensorAttribute(std::string s, Mat3d*  ta);
		void setTensorAttributeZero(std::string s);
		void setTensorAttributeIdentity(std::string s);
		void setDpDensity(double dp, double rho);
		void setDpActual(double dp) { dp_ = dp; }
		void create();
		double getDp() const { return dp_; }
	};
}

#endif // !_ENSEMBLE_CREATOR_