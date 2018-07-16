/*! \file m_ensemble_creator.cpp */

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
#include "m_ensemble_creator.h"
#include <iostream>

namespace msl {

	//==============================================================================
	/// creator routines: 
	/// creator constructor-> addAttributes-> setOrientation/Dims
	/// -> setDp -> create -> setAttribute
	//==============================================================================

	//==============================================================================
	/// Ctor
	//==============================================================================
	EnsembleCreator::EnsembleCreator(std::string name, Vec3d offset) 
		: cache_(0), offset_(offset)
	{
		ensemble_ = std::make_shared<Ensemble>();
		shape_ = ShapeFactory::ShapeBuilder(name);
		shape_actual_ = ShapeFactory::ShapeBuilder(name);
		class_name_ = "EnsembleCreator";
	}

	//==============================================================================
	/// Destructor
	//==============================================================================
	EnsembleCreator::~EnsembleCreator()	{
		if (cache_ != 0)
			delete[] cache_;
	}

	//==============================================================================
	/// Adding scalar attributes
	//==============================================================================
	void EnsembleCreator::addScalarAttributes(std::vector<std::string> scalars)	{
		for (auto sa : scalars)
			ensemble_->addScalarAttribute(sa);
	}

	//==============================================================================
	/// Adding vector attributes
	//==============================================================================
	void EnsembleCreator::addVectorAttributes(std::vector<std::string> vectors) {
		for (auto va : vectors)
			ensemble_->addVectorAttribute(va);
	}

	//==============================================================================
	/// Adding tensor attribute
	//==============================================================================
	void EnsembleCreator::addTensorAttributes(std::vector<std::string> tensors) {
		for (auto ta : tensors)
			ensemble_->addTensorAttribute(ta);
	}

	//==============================================================================
	/// Set scalar attribute values
	//==============================================================================
	void EnsembleCreator::setScalarAttribute(std::string s, double * sa) {
		int np = ensemble_->getNp();
		Ensemble::ScalarAttrPtr sa_ptr = ensemble_->getScalarAttrPtr(s);
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np; i++)
			sa_ptr->getAttr()[i] = sa[i];
	}

	//==============================================================================
	/// Set scalar attribute as constant
	//==============================================================================
	void EnsembleCreator::setScalarAttributeConstant(std::string s, double c) {
		int np = ensemble_->getNp();
		Ensemble::ScalarAttrPtr sa_ptr = ensemble_->getScalarAttrPtr(s);
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np; i++)
			sa_ptr->getAttr()[i] = c;
	}

	//==============================================================================
	/// Set vector attribute
	//==============================================================================
	void EnsembleCreator::setVectorAttribute(std::string s, Vec3d * va) {
		int np = ensemble_->getNp();
		Ensemble::VectorAttrPtr va_ptr = ensemble_->getVectorAttrPtr(s);
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np; i++)
			va_ptr->getAttr()[i] = va[i];
	}

	//==============================================================================
	/// Set vector attribute as zeros
	//==============================================================================
	void EnsembleCreator::setVectorAttributeZero(std::string s) {
		int np = ensemble_->getNp();
		Ensemble::VectorAttrPtr va_ptr = ensemble_->getVectorAttrPtr(s);
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np; i++)
			va_ptr->getAttr()[i] = Vec3d::Zero();
	}

	//==============================================================================
	/// Set tensor attribute
	//==============================================================================
	void EnsembleCreator::setTensorAttribute(std::string s, Mat3d * ta) {
		int np = ensemble_->getNp();
		Ensemble::TensorAttrPtr ta_ptr = ensemble_->getTensorAttrPtr(s);
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np; i++)
			ta_ptr->getAttr()[i] = ta[i];
	}

	//==============================================================================
	/// Set tensor attribute as zeros
	//==============================================================================
	void EnsembleCreator::setTensorAttributeZero(std::string s) {
		int np = ensemble_->getNp();
		Ensemble::TensorAttrPtr ta_ptr = ensemble_->getTensorAttrPtr(s);
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np; i++)
			ta_ptr->getAttr()[i] = Mat3d::Zero();
	}

	//==============================================================================
	/// Set tensor attribute as identity matrix
	//==============================================================================
	void EnsembleCreator::setTensorAttributeIdentity(std::string s) {
		int np = ensemble_->getNp();
		Ensemble::TensorAttrPtr ta_ptr = ensemble_->getTensorAttrPtr(s);
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np; i++)
			ta_ptr->getAttr()[i] = Mat3d::Identity();
	}

	//==============================================================================
	/// set dp
	//==============================================================================
	void EnsembleCreator::setDpDensity(double dp, double rho) {
		if (shape_->getVolume() == 0.0)
			throwException("create()", "Must set dims first");
		double l = shape_->getFirstDim();
		int sz = int(l / dp) + 1;
		dp_ = l / double(sz);
		density_ = rho;
	}

	//==============================================================================
	/// create ensemble
	//==============================================================================
	void EnsembleCreator::create() {
		if (dp_ == 0.0)
			throwException("create()", "Must set dp first");
		double voxel = shape_->shape2D() ? dp_ * dp_ : dp_ * dp_ * dp_;
		std::cout << "voxel volume: " << voxel << std::endl;
		std::cout << "shape volume: " << shape_->getVolume() << std::endl;
		int size_approx = int(shape_->getVolume() / voxel * 1.2);
		std::cout << "size approx: " << size_approx << std::endl;
		if (cache_ != 0)
			delete[] cache_;
		try {
			cache_ = new Vec3d[size_approx];
		}
		catch (std::bad_alloc) {
			throwException("create()", "Error occured allocating memory for cache");
		}
		std::vector<double> dims = shape_->getDims();
		Vec3i nums(0, 0, 0);
		for (int i = 0; i < 3; i++)
			nums[i] = round(dims[i] / dp_);
		int cnt = 0;
		Vec3d origin = dp_ * 0.5 * (nums - Vec3i::Ones()).cast<double>();
		Mat3d rotation = shape_->getRotation();
		if (shape_->shape2D()) { 
			origin[1] = 0.0; 
			for (int i = 0; i < nums[0]; i++) {
				for (int j = 0; j < nums[2]; j++) {
					Vec3d pos = Vec3d(double(i) * dp_, 0.0, double(j) * dp_) - origin;
					if (shape_->isWithin(pos))
						cache_[cnt++] = rotation * pos + offset_;
				}
			}
		}
		else {
			for (int i = 0; i < nums[0]; i++) {
				for (int j = 0; j < nums[1]; j++) {
					for (int k = 0; k < nums[2]; k++) {
						Vec3d pos = dp_ * Vec3d(double(i), double(j), double(k)) - origin;
						if (shape_->isWithin(pos))
							cache_[cnt++] = rotation * pos + offset_;
					}
				}
			}		 
		}
		ensemble_->setNp(cnt);
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < cnt; i++) {
			(ensemble_->pos_)[i] = cache_[i];
			(ensemble_->acc_)[i] = Vec3d::Zero();
			(ensemble_->vel_)[i] = Vec3d::Zero();
		}
		ensemble_->setDpDensity(dp_, density_);
		if (ensemble_->hasVectorAttribute("Initial_position"))
			setVectorAttribute("Initial_position", cache_);
	}
}
