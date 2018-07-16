/*! \file m_ensemble.h */

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
#ifndef _M4_ENSEMBLE_
#define _M4_ENSEMBLE_

#include "m_ensemble_attribute.h"
#include <vector>
#include <memory>
#include <unordered_map>
#pragma warning (disable:4244)
#pragma warning (disable:4267)

namespace msl {
	class Ensemble : protected MsObj {
	public:
		friend class EnsembleCreator;
		typedef std::shared_ptr<ScalarAttribute> ScalarAttrPtr;
		typedef std::shared_ptr<VectorAttribute> VectorAttrPtr;
		typedef std::shared_ptr<TensorAttribute> TensorAttrPtr;
		typedef std::unordered_map<std::string, int> Code;

		Ensemble();

		Ensemble(const Ensemble& en) = delete;
		
		Ensemble& operator=(const Ensemble& en) = delete;

		~Ensemble();

		void setNp(int i);
		void setDpDensity(double dp, double rho);
		void addScalarAttribute(std::string s = "");
		void addVectorAttribute(std::string s = "");
		void addTensorAttribute(std::string s = "");
		int getNp() const { return np_; }
		Vec3d getPosMin() const { return pos_min_; }
		Vec3d getPosMax() const { return pos_max_; }
		double getDp() const { return dp_; }
		double getDensity() const { return density_; }
		double getMass() const { return mass_; }
		Vec2d getVelMinMax() const;
		Vec2d getAccMinMax() const;
		int* getId() { return id_; }
		int* getDict() { return dict_; }
		Vec3d* getPos() { return pos_; }
		Vec3d* getVel() { return vel_; }
		Vec3d* getAcc() { return acc_; }
		void calcBox();

		bool hasScalarAttribute(std::string s) { return scalar_.find(s) != scalar_.end(); }
		bool hasVectorAttribute(std::string s) { return vector_.find(s) != vector_.end(); }
		bool hasTensorAttribute(std::string s) { return tensor_.find(s) != tensor_.end(); }
		ScalarAttrPtr getScalarAttrPtr(std::string s);
		VectorAttrPtr getVectorAttrPtr(std::string s);
		TensorAttrPtr getTensorAttrPtr(std::string s);
		std::vector<ScalarAttrPtr> getScalarAttrPtrs() { return scalar_attrs_; }
		std::vector<VectorAttrPtr> getVectorAttrPtrs() { return vector_attrs_; }
		std::vector<TensorAttrPtr> getTensorAttrPtrs() { return tensor_attrs_; }
		double getMemoryInMB() const;

	private:
		int			np_;		//! number of total particles
		int*		id_;		//! particle ids
		int*		dict_;		//! used to sort particles  pid-> curr
		double		density_;	//! density
		double		mass_;		//! mass of single particle
		double		dp_;		//! inter-particle distance
		Vec3d*		pos_;		//! position vectors for ensemble
		Vec3d*		vel_;		//! velocity vectors for ensemble
		Vec3d*		acc_;		//! acceleration vectors for ensemble
		Vec3d		pos_min_;	//! box containing all particles pos min
		Vec3d		pos_max_;	//! box containing all particles pos max
		Code		scalar_;	//! map to scalar attributes
		Code		vector_;	//! map to vector attributes
		Code		tensor_;	//! map to tensor attributes
		std::vector<ScalarAttrPtr>	scalar_attrs_;
		std::vector<VectorAttrPtr>	vector_attrs_;
		std::vector<TensorAttrPtr>  tensor_attrs_;
	};
}

#endif // ! _M4_ENSEMBLE_