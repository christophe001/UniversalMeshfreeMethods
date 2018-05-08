/*! \file m_ensemble.cpp */

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

#include "m_ensemble.h"
#include <algorithm>

namespace msl {

	//==============================================================================
	/// Ctor
	//==============================================================================
	Ensemble::Ensemble() :
		np_(0), id_(0), dict_(0),
		pos_(0), vel_(0), acc_(0)
	{
		class_name_ = "Ensemble";
		pos_min_ = Vec3d::Constant(std::numeric_limits<double>::max());
		pos_max_ = Vec3d::Constant(std::numeric_limits<double>::min());
	}

	//==============================================================================
	/// Destructor
	//==============================================================================
	Ensemble::~Ensemble() {
		if (id_ != 0)			delete[] id_;
		if (dict_ != 0)			delete[] dict_;
		if (pos_ != 0)			delete[] pos_;
		if (vel_ != 0)			delete[] vel_;
		if (acc_ != 0)			delete[] acc_;
	}

	//==============================================================================
	/// Set memory for np particles
	//==============================================================================
	void Ensemble::setNp(int np) {
		if (np <= 0)
			throwException("setNp", "Invalid particle numbers");
		if (id_ != 0)		delete[] id_;
		if (dict_ != 0)		delete[] dict_;
		if (pos_ != 0)		delete[] pos_;
		if (vel_ != 0)		delete[] vel_;
		if (acc_ != 0)		delete[] acc_;

		np_ = np;
		try {
			id_ = new int[np];
			dict_ = new int[np];
			for (int i = 0; i < np; i++) { id_[i] = i; dict_[i] = i; }
			pos_ = new Vec3d[np];
			vel_ = new Vec3d[np];
			acc_ = new Vec3d[np];
		}
		catch (std::bad_alloc) {
			throwException("setNp", "Error occured while allocating memory");
		}

		for (auto&& sa : scalar_attrs_) sa->setNp(np);
		for (auto&& va : vector_attrs_) va->setNp(np);
		for (auto&& ta : tensor_attrs_) ta->setNp(np);
	}

	//==============================================================================
	/// Add scalar attribute to ensemble
	//==============================================================================
	void Ensemble::addScalarAttribute(std::string s) {
		int sz = scalar_attrs_.size();
		auto sa = std::make_shared<ScalarAttribute>(s);
		scalar_attrs_.push_back(std::move(sa));
		scalar_[s] = sz;
	}

	//==============================================================================
	/// Add vector attribute to ensemble
	//==============================================================================
	void Ensemble::addVectorAttribute(std::string s) {
		int sz = vector_attrs_.size();
		auto va = std::make_shared<VectorAttribute>(s);
		vector_attrs_.push_back(std::move(va));
		vector_[s] = sz;
	}

	//==============================================================================
	/// Add tensor(matrix) attribute to ensemble
	//==============================================================================
	void Ensemble::addTensorAttribute(std::string s) {
		int sz = tensor_attrs_.size();
		auto ta = std::make_shared<TensorAttribute>(s);
		tensor_attrs_.push_back(std::move(ta));
		tensor_[s] = sz;
	}

	//==============================================================================
	/// Get max vel/ min vel magnitude
	//==============================================================================
	Vec2d Ensemble::getVelMinMax() const {
		Vec2d res(std::numeric_limits<double>::max(), 0);
		for (int i = 0; i < np_; i++) {
			res[0] = std::min(vel_[i].norm(), res[0]);
			res[1] = std::max(vel_[i].norm(), res[1]);
		}
		return res;
	}

	//==============================================================================
	/// get max acc/ min acc magnitude
	//==============================================================================
	Vec2d Ensemble::getAccMinMax() const {
		Vec2d res(std::numeric_limits<double>::max(), 0);
		for (int i = 0; i < np_; i++) {
			res[0] = std::min(acc_[i].norm(), res[0]);
			res[1] = std::max(acc_[i].norm(), res[1]);
		}
		return res;
	}

	//==============================================================================
	/// Calculate minimized axis aligned box containing ensemble
	//==============================================================================
	void Ensemble::calcBox() {
		for (int i = 0; i < np_; i++) {
			pos_min_ = pos_min_.cwiseMin(pos_[i]);
			pos_max_ = pos_max_.cwiseMax(pos_[i]);
		}
	}

	//==============================================================================
	/// Retrieve pointer to desired scalar attribute
	//==============================================================================
	Ensemble::ScalarAttrPtr Ensemble::getScalarAttrPtr(std::string s) {
		if (scalar_.find(s) == scalar_.end())
			throwException("getScalarAttrPtrs", "Cannot find attribute: " + s);
		return scalar_attrs_[scalar_[s]];
	}

	//==============================================================================
	/// Retrieve pointer to desired vector attribute
	//==============================================================================
	Ensemble::VectorAttrPtr Ensemble::getVectorAttrPtr(std::string s) {
		if (vector_.find(s) == vector_.end())
			throwException("getVectorAttrPtrs", "Cannot find attribute: " + s);
		return vector_attrs_[vector_[s]];
	}

	//==============================================================================
	/// Retrieve pointer to desired tensor(matrix) attribute
	//==============================================================================
	Ensemble::TensorAttrPtr Ensemble::getTensorAttrPtr(std::string s) {
		if (tensor_.find(s) == tensor_.end())
			throwException("getTensorAttrPtrs", "Cannot find attribute: " + s);
		return tensor_attrs_[tensor_[s]];
	}

	//==============================================================================
	/// Get memory used for Ensemble in MB
	//==============================================================================
	double Ensemble::getMemoryInMB() const {
		int64_t size = sizeof(int) * (2 * np_ + 1) + sizeof(Vec3d) * np_ * 3
			+ sizeof(Vec3d*) * 3 + sizeof(int*) * 2 + sizeof(std::vector<ScalarAttrPtr>)
			+ sizeof(std::vector<VectorAttrPtr>) + sizeof(std::vector<TensorAttrPtr>);
		double memsize = size / 1048576.0;
		for (auto&& sa : scalar_attrs_)
			memsize += sa->getMemoryInMB();
		for (auto&& va : vector_attrs_)
			memsize += va->getMemoryInMB();
		for (auto&& ta : tensor_attrs_)
			memsize += ta->getMemoryInMB();
		return memsize;
	}

}