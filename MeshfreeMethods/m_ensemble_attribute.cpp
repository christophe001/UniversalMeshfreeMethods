/*! \file m_ensemble_attribute.cpp */

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

#include "m_ensemble_attribute.h"

namespace msl {
	//******************************************************************************
	/// class BondAttribute
	//******************************************************************************
	//==============================================================================
	/// Set bond list size(approximate)
	//==============================================================================
	void BondAttribute::setListSize(long ls) {
		if (ls <= 0) throwException("setListSize", "Invalid list size");
		if (attr_ != 0) delete[] attr_;
		list_size_ = ls;
		try {
			attr_ = new double[list_size_];
		}
		catch (std::bad_alloc) {
			throwException("setListSize", "Error occured while allocating memeory");
		}
	}

	//==============================================================================
	///	Get memory used for BondAttribute
	//==============================================================================
	double BondAttribute::getMemoryInMB() const
	{
		int64_t size = sizeof(long) + sizeof(double) * list_size_ + sizeof(double*);
		double memsize = size / 1048576.0;
		return memsize;
	}

	//******************************************************************************
	/// class ScalarAttribute
	//******************************************************************************
	//==============================================================================
	/// Set scalar attribute size
	//==============================================================================
	void ScalarAttribute::setNp(int np) {
		if (np <= 0) throwException("setNp", "Invalid particle number");
		if (attr_ != 0) delete[] attr_;
		np_ = np;
		try {
			attr_ = new double[np];
		}
		catch (std::bad_alloc) {
			throwException("setNp", "Error occured while allocating memory");
		}
	}

	//==============================================================================
	///	Get memory used for ScalarAttribute
	//==============================================================================
	double ScalarAttribute::getMemoryInMB() const {
		int64_t size = sizeof(int) + sizeof(double) * np_ + sizeof(double*);
		double memsize = size / 1048576.0;
		return memsize;
	}


	//******************************************************************************
	/// class VectorAttribute
	//******************************************************************************
	//==============================================================================
	/// Set vector attribute size
	//==============================================================================
	void VectorAttribute::setNp(int np) {
		if (np <= 0) throwException("setNp", "Invalid particle number");
		if (attr_ != 0) delete[] attr_;
		np_ = np;
		try {
			attr_ = new Vec3d[np];
		}
		catch (std::bad_alloc) {
			throwException("setNp", "Error occured while allocating memory");
		}
	}

	//==============================================================================
	///	Get memory used for VectorAttribute
	//==============================================================================
	double VectorAttribute::getMemoryInMB() const {
		int64_t size = sizeof(int) + sizeof(Vec3d) * np_ + sizeof(Vec3d*);
		double memsize = size / 1048576.0;
		return memsize;
	}

	//******************************************************************************
	/// class TensorAttribute
	//******************************************************************************
	//==============================================================================
	/// Set tensor attribute size
	//==============================================================================
	void TensorAttribute::setNp(int np) {
		if (np <= 0) throwException("setNp", "Invalid particle number");
		if (attr_ != 0) delete[] attr_;
		np_ = np;
		try {
			attr_ = new Mat3d[np];
		}
		catch (std::bad_alloc) {
			throwException("setNp", "Error occured while allocating memory");
		}
	}

	//==============================================================================
	///	Get memory used for TensorAttribute
	//==============================================================================
	double TensorAttribute::getMemoryInMB() const {
		int64_t size = sizeof(int) + sizeof(Mat3d) * np_ + sizeof(Mat3d*);
		double memsize = size / 1048576.0;
		return memsize;
	}
}