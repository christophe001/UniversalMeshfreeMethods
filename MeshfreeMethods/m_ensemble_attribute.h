/*! \file m_ensemble_attribute.h*/

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

#ifndef _M4_ENSEMBLE_ATTR_
#define _M4_ENSEMBLE_ATTR_

#include "m_types.h"
#include "m_obj.h"
#include <string>
#pragma warning (disable:4244)
#pragma warning (disable:4267)

namespace msl {

	class BondAttribute : protected MsObj {
	private:
		long		list_size_;		//! Bond list size, length of attribute
		double*		attr_;			//! Pointer to attribute
		std::string info_;			//! Infomation of attribute

	public:

		//! Ctor
		BondAttribute(std::string s = "") : list_size_(0), attr_(0) {
			class_name_ = "BondAttribute";
			info_ = s;
		}

		//! Destructor
		~BondAttribute() { if (attr_ != 0) delete[] attr_; }

		long getListSize() const { return list_size_; }

		//! Return attribute pointer
		double* getAttr() { return attr_; }

		//! Return list size
		void setListSize(long ls);

		//! Set infomation for attribute
		void setInfo(std::string s) { info_ = s; }

		//! Return info
		std::string getInfo() const { return info_; }

		double getMemoryInMB() const;
	};

	class ScalarAttribute : protected MsObj {
	private:
		int			np_;		//! Num of points, length of attribute
		double*		attr_;		//! Pointer to attribute
		std::string	info_;		//! Infomation of attribute

	public:
		//! Ctor
		ScalarAttribute(std::string s = "") : np_(0), attr_(0) {
			class_name_ = "ScalarAttribute";
			info_ = s;
		}

		//! Destructor
		~ScalarAttribute() { if (attr_ != 0) delete[] attr_; }

		//! Return size
		int getNp() const { return np_; }

		//! Return pointer to attribute
		double*  getAttr() { return attr_; }

		//! Set size
		void setNp(int np);

		// !Set infomation for attribute
		void setInfo(std::string s) { info_ = s; }

		//! Return info
		std::string getInfo() const { return info_; }

		double getMemoryInMB() const;
	};

	class VectorAttribute : protected MsObj {
	private:
		int			np_;
		Vec3d*		attr_;
		std::string	info_;

	public:
		VectorAttribute(std::string s = "") : np_(0), attr_(0) {
			class_name_ = "VectorAttribute";
			info_ = s;
		}
		~VectorAttribute() { if (attr_ != 0) delete[] attr_; }

		int getNp() const { return np_; }

		Vec3d*  getAttr() { return attr_; }

		void setNp(int N);

		void setInfo(std::string s) { info_ = s; }

		std::string getInfo() const { return info_; }

		double getMemoryInMB() const;
	};

	class TensorAttribute : protected MsObj {
	private:
		int			np_;
		Mat3d*		attr_;
		std::string	info_;

	public:
		TensorAttribute(std::string s = "") : np_(0), attr_(0) {
			class_name_ = "TensorAttribute";
			info_ = s;
		}

		~TensorAttribute() { if (attr_ != 0) delete[] attr_; }

		int getNp() const { return np_; }

		Mat3d*  getAttr() { return attr_; }

		void setNp(int N);

		void setInfo(std::string s) { info_ = s; }

		std::string getInfo() const { return info_; }

		double getMemoryInMB() const;
	};
}

#endif // !_M4_ENSEMBLE_ATTR_
