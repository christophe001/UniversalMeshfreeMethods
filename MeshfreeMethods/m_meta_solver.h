/*! \file m_meta_solver.h */

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

#ifndef _M4_META_SOLVER_
#define _M4_META_SOLVER_

#include "m_meshfree_base.h"
#include "m_contact_manager.h"
#include "m_contact_manager_rigid_body.h"

namespace msl {
	
	class MetaSolver : public MsObj {
	protected:
		std::shared_ptr<MeshfreeBase>		slave_solv_;
		std::shared_ptr<SortEnsemble>		slave_sort_;
		std::shared_ptr<MeshfreeBase>		master_solv_;
		std::shared_ptr<SortEnsemble>		master_sort_;
		std::shared_ptr<ContactMngerSRB>	srb_cm_;
		bool								use_srb_;
		DomainConfig						domain_config_;
		std::shared_ptr<ContactManager>		cm_;
		double								dp_s_, dp_m_;
		ObjInfo								master_obj_, slave_obj_;

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
		//! relaxation parameter
		double									rc_;	//! relaxation parameter

		//*********************************************************************
		//							I/O parameters
		//*********************************************************************
		//! ensemble attributes to save
		std::vector<std::string>				s_scalar_saves_;
		std::vector<std::string>				m_scalar_saves_;
		std::vector<std::string>				s_vector_saves_;
		std::vector<std::string>				m_vector_saves_;
		//! file name and target folder
		std::string								filename_;
		std::string								folder_;
		std::string								slave_name_, master_name_;

	public:
		//*********************************************************************
		//							Running routine
		//*********************************************************************
		// 1. MetaSolver(const_timesetp, enable_relax, rc_param)
		// 2. configIO()
		// 3. configRun()
		// 4. createScene()
		// 5. configContact()
		// 6. conditions
		///   setInitVel()
		///   configPosConstraint(shape)				
		///   configExternalPosForce(shape, vec)
		///   configExternalInitPosForce(shape, vec)
		// 7. run()
		
		MetaSolver(bool ct = true, bool er = false, double rc = 0.0);
		virtual ~MetaSolver() {}
		void setDomainConfig(const DomainConfig& domain_config) { domain_config_ = domain_config; }
		void printLabel() const;
		void printInfo() const;
		void printMemoryInfo() const;
		virtual void configMeshfreeSolver(std::string slav_lc, std::string slav_ec, std::vector<std::string> slav_attrs,
		std::string master_lc, std::string master_ec, std::vector<std::string> master_attrs);
		virtual void configContact();
		void configSRB();
		void createScene(ObjInfo master, ObjInfo slave, DomainConfig domain_config);
		void configRun(double dt, double total_time, int sv_step);
		void addConfinedRegions(std::vector<std::shared_ptr<Shape>> regions);
		void addConfinedRegionsComplement(std::vector<std::shared_ptr<Shape>> regions);
		void run();
		void runWithNoSlip();
		void saveVtk();
		void saveVtkAll();
		void saveVtkSlave();
		void saveVtkMaster();
		void verletUpdate();
		void verletUpdateWithNoSlip();
		void configIO(std::string folder, std::string filename,
			std::string slave = "target", std::string master = "projectile");
	};
}

#endif // !_M4_META_SOLVER_

