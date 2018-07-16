/*! \file m_meshfree_solver.h */

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

#ifndef _M4_MESHFREE_
#define _M4_MESHFREE_

#include <iostream>
#include "m_lagrangian_compute.h"
#include "m_state_kl_model.h"
#include "m_rigid_body.h"
#include "m_eulerian_compute.h"
#include "m_shape.h"
#include "m_ensemble_creator.h"
//#include "m_state_pd_jh2.h"

namespace msl {
	struct ObjInfo {
		std::string shape;
		Vec3d		dims;
		Vec3d		offset;
		Vec3d       orientation;
		double		theta;
		double		dp;
		double      density;
		Vec3d	    velocity;
		std::string solver;
		ObjInfo(std::string s, Vec3d d, std::string solv, double delta_d, double rho, Vec3d vel = Vec3d::Zero(),
			Vec3d os = Vec3d::Zero(), Vec3d or = Vec3d::UnitZ(), double th = 0.0) {
			shape = s;
			dims = d;
			solver = solv;
			dp = delta_d;
			density = rho;
			velocity = vel;
			offset = os;
			orientation = or ;
			theta = th;
		}
		ObjInfo& operator=(const ObjInfo& obj) {
			shape = obj.shape;
			dims = obj.dims;
			offset = obj.offset;
			dp = obj.dp;
			density = obj.density;
			velocity = obj.velocity;
			orientation = obj.orientation;
			theta = obj.theta;
			solver = obj.solver;
			return *this;
		}
	};

	struct AttributeInfo {
		std::string m_attr_type;
		std::string m_name;
		bool		m_save;
		std::string m_init;
		double		m_val;
		AttributeInfo() : m_attr_type("scalar"), m_save(false), m_name("_"), m_init("zero"), m_val(0) {}
		AttributeInfo(std::string attr, std::string name, bool save, std::string init, double val = 0)
			: m_attr_type(attr), m_name(name), m_save(save), m_init(init), m_val(val) {}
		AttributeInfo& operator=(const AttributeInfo& info) {
			m_attr_type = info.m_attr_type;
			m_name = info.m_name;
			m_save = info.m_save;
			m_init = info.m_init;
			m_val = info.m_val;
			return *this;
		}
	};

	class MeshfreeBase : public MsObj {
	protected:
		//*********************************************************************
		//		Ensemble structure and computational features
		//*********************************************************************
		//! ensemble creator
		std::shared_ptr<EnsembleCreator>		creator_;
		//! underlying ensemble
		std::shared_ptr<Ensemble>				ensemble_;
		//! pointer to sorted ensemble
		std::shared_ptr<SortEnsemble>			sorted_;
		//! neighborhood data
		std::shared_ptr<NeighborhoodData>		nbh_;
		//! peri-neighborhood data
		std::shared_ptr<PeriNeighborData>		pbh_;
		//! compute object manage bond based compute
		std::shared_ptr<LagrangianCompute>		l_compute_;
		//! compute object manage eulerian compute
		std::shared_ptr<EulerianCompute>		e_compute_;
		//! domain config
		DomainConfig							domain_config_;
		
		//*********************************************************************
		//				Ensemble properties and attributes
		//*********************************************************************
		//! number of total particles
		int										np_;
		int*									id_;
		int*									dict_;
		//! inter-particle distance
		double									dp_;
		//! material density
		double									density_;
		//! horizon for particle interaction
		double									horizon_;
		//! key attributes, pos, vel and acc
		Vec3d*									pos_;
		Vec3d*									vel_;
		Vec3d*									acc_;
		//! memory stack for storage
		Vec3d*									stack_;
		//! ensemble attributes
		std::vector<Ensemble::ScalarAttrPtr>	scalar_attrs_;
		std::vector<Ensemble::VectorAttrPtr>	vector_attrs_;
		std::vector<Ensemble::TensorAttrPtr>	tensor_attrs_;
		//! ensemble bond attribute names
		std::vector<std::string>				bond_attrs_;

		//*********************************************************************
		//							I/O parameters
		//*********************************************************************
		//! ensemble attributes to save
		std::vector<std::string>				scalar_saves_;
		std::vector<std::string>				vector_saves_;
		//! file name and target folder
		std::string								filename_;
		std::string								folder_;

		//*********************************************************************
		//!							Run parameters
		//*********************************************************************
		//! use of constant timestep
		bool									constant_timestep_;
		//! current timestep
		int										timestep_;
		//! save info every sv_step_ timestep (constant ts), current save count 
		int										sv_step_, sv_count_;
		//! timestep when setting constant_timestep_ true
		double									dt_max_;
		//! actual timestep, total simulation time, actual sim time
		double									dt_, T_, sim_time_;
		//! sim time since last save
		double									t_cur_;
		//! time duration between saves
		double									t_intv_;
		//! if adopt adaptive relaxation
		bool									relaxation_;
		//! if adopt viscosity
		bool									viscosity_;
		//! relaxation parameter
		double									rc_;	//! relaxation parameter

		//*********************************************************************
		//			Force and velocity initialization parameters
		//*********************************************************************
		//! forces applied based on current position
		std::vector<std::shared_ptr<Shape>>		force_regions_;
		std::vector<Vec3d>						forces_;
		//! confined region for position constraint
		std::vector<std::shared_ptr<Shape>>		confined_regions_;
		//! forces applied based on initial position
		std::vector<std::shared_ptr<Shape>>		force_init_regions_;
		std::vector<Vec3d>						forces_init_;
		//! initial velocity
		Vec3d									init_vel_;

	public:
		MeshfreeBase(bool const_ts = true, bool enable_rx = false, double rc = 0.0);

		void createEnsemble(ObjInfo obj, DomainConfig domain_config, std::vector<AttributeInfo>& infos);
		virtual ~MeshfreeBase() {}
		//void setParams(double horizon, double sv_step, bool relaxation, double rc = 0.0);
		void setDomainConfig(const DomainConfig& domain_config) { domain_config_ = domain_config; }
		void printLabel() const;
		void printInfo() const;
		void printMemoryInfo() const;
		void configSolver(std::string lcompute, std::string ecompute,
			std::vector<std::string> attrs);
		//! config state kl model
		void configLagrangianKL();
		//! config rigid body model
		void configRigidBody();
		void adaptiveRelaxation();
		void setInitVel(Vec3d vel);
		void addPosForce(std::shared_ptr<Shape> shape, Vec3d force);
		void addPosConstraint(std::shared_ptr<Shape> shape);
		void addInitialPosForce(std::shared_ptr<Shape> shape, Vec3d force);
		void configExternalPosForce(std::shared_ptr<Shape> shape, Vec3d force);
		void configExternalInitPosForce(std::shared_ptr<Shape> shape, Vec3d force);
		void initForces();
		void setDt(double dt);
		double getDt() const { return dt_; }
		double getSaveDt() const { return t_intv_; }
		virtual void configDt() {}
		double getDtMax() const { return dt_max_; }	
		void addViscosity() {}
		void verletUpdate();
		void saveVtk();
		void configIO(std::string folder, std::string filename);
		void configRun(double dt, double total_time, int sv_step);
		virtual void run();
		friend class MetaSolver;
	};
}
#endif
