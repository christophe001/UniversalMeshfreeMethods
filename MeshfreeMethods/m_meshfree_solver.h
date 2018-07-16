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

#ifndef _M4_MESHFREE_SOLVER_
#define _M4_MESHFREE_SOLVER_

#include <iostream>
#include "m_lagrangian_compute.h"
#include "m_eulerian_compute.h"
#include "m_state_pd_jh2.h"

namespace msl {
	class MeshfreeSolver : public MsObj {
	protected:
		//typedef std::vector<std::shared_ptr<LagrangianCompute>> LCompute;
		//typedef std::vector<std::shared_ptr<EulerianCompute>>	ECompute;

		//*********************************************************************
		//!							Underlying ensemble
		//*********************************************************************
		std::shared_ptr<Ensemble>				ensemble_;
		std::shared_ptr<SortEnsemble>			sorted_;
		std::shared_ptr<NeighborhoodData>		nbh_;
		std::shared_ptr<PeriNeighborData>		pbh_;
		
		//*********************************************************************
		//!				Lagrangian and Eulerian compute object 
		//*********************************************************************
		std::shared_ptr<LagrangianCompute>		l_compute_;
		std::shared_ptr<EulerianCompute>		e_compute_;

		//*********************************************************************
		//!						creator of ensemble
		//*********************************************************************
		std::shared_ptr<EnsembleCreator>		creator_;

		//*********************************************************************
		//!							ensemble attributes
		//*********************************************************************
		DomainConfig							domain_config_;
		int										np_;
		int*									id_;
		int*									dict_;
		Vec3d*									pos_;
		Vec3d*									vel_;
		Vec3d*									acc_;
		//! Memory used to update ensemble attributes
		Vec3d*									stack_; 

		//*********************************************************************
		//!			Ensemble attributes associated with particles
		//*********************************************************************
		std::vector<Ensemble::ScalarAttrPtr>	scalar_attrs_;
		std::vector<Ensemble::VectorAttrPtr>	vector_attrs_;
		std::vector<Ensemble::TensorAttrPtr>	tensor_attrs_;

		//*********************************************************************
		//!				Attributes to be saved during simulation
		//*********************************************************************
		std::vector<std::string>				scalar_saves_;
		std::vector<std::string>				vector_saves_;

		//*********************************************************************
		//!							Run parameters
		//*********************************************************************
		std::string								filename_;
		std::string								folder_;
		int										timestep_, max_step_;
		int										sv_step_, sv_count_;
		double									dp_, dt_, T_;
		double									dt_max_;
		double									horizon_;
		bool									relaxation_;
		double									rc_;	//! relaxation parameter

		//*********************************************************************
		//!							Initial forces
		//*********************************************************************
		std::vector<std::shared_ptr<Shape>>		force_regions_;
		std::vector<std::shared_ptr<Shape>>		force_init_regions_;
		std::vector<Vec3d>						forces_;
		std::vector<Vec3d>						forces_init_;
		

	public:
		MeshfreeSolver();
		class AttributeInfo {
		public:
			std::string m_attr_type;
			std::string m_name;
			bool		m_save;
			std::string m_init;
			double		m_val;
			AttributeInfo() : m_attr_type("scalar"), m_save(false), m_name("_"), m_init("zero"), m_val(0) {}
			AttributeInfo(std::string attr, std::string name, bool save, std::string init, double val = 0)
				: m_attr_type(attr), m_name(name), m_save(save), m_init(init), m_val(val){}
			AttributeInfo& operator=(const AttributeInfo& info) {
				m_attr_type = info.m_attr_type;
				m_name = info.m_name;
				m_save = info.m_save;
				m_init = info.m_init;
				m_val = info.m_val;
				return *this;
			}
		};
		std::vector<AttributeInfo> getJH2AttributeConfig(double lambda, double mu);

		virtual ~MeshfreeSolver();
		void setParams(double horizon, double sv_step, bool relaxation, double rc = 0.0);
		void setJH2Params(int i);
		void setDomainConfig(DomainConfig domain_config) { domain_config_ = domain_config; }
		void printLabel();
		void printInfo();
		void createEnsemble(std::string shape, DomainConfig domain_config, std::vector<AttributeInfo>& infos, 
							Vec3d dims, Vec3d offset, Vec3d orientation , double theta, double dp);
		void configSolver(std::string lcompute, std::string ecompute, 
			std::vector<std::string> attrs = std::vector<std::string>());
		void adaptiveRelaxation();
		void setInitVel(Vec3d vel);
		void addPosForce(std::shared_ptr<Shape> shape, Vec3d force);
		void addInitialPosForce(std::shared_ptr<Shape> shape, Vec3d force);
		void configExternalPosForce(std::shared_ptr<Shape> shape, Vec3d force);
		void configExternalInitPosForce(std::shared_ptr<Shape> shape, Vec3d force);
		void initForces();
		void setDt(double dt) { dt_ = dt; }
		void printMemoryInfo();
		void verletUpdate();
		void saveVtk();
		void configRunParams(std::string folder, std::string filename, double dt, double total_time);
		void run();
	};
}

#endif