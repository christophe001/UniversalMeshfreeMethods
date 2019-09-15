/*! \file m_meta_solver.cpp */

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

#include "m_meta_solver.h"
#include "u_timer.h"
#include <windows.h>
#include <fstream>
#include <cmath>

constexpr auto theWatch = msl::Timer::instance;

namespace msl {
	MetaSolver::MetaSolver(bool ct, bool er, double rc) {
		class_name_ = "MetaSolver";
		slave_solv_ = std::make_shared<MeshfreeBase>(3.5);
		master_solv_ = std::make_shared<MeshfreeBase>(3.5);
		constant_timestep_ = ct;
		relaxation_ = er;
		rc_ = rc;
		use_srb_ = false;
	}

	void MetaSolver::printLabel() const {
		std::cout
			<< "************************************************************************\n"
			<< "\n"
			<< "         ___    ___      _ _   _               _          ___			\n"
			<< "         |  \\  /  |     | | | (_)             | |        /   |			\n"
			<< "         | . \\/ . |_   _| | |_ _ ___  ___ __ _| | ___   / /| |			\n"
			<< "         | |\\  /| | | | | | __| / __|/ __/ _` | |/ _ \\ / /_| |	    \n"
			<< "         | | \\/ | | |_| | | |_| \\__ \\ (_| (_| | |  __/ \\___  |		\n"
			<< "         \\_|    |_/\\__,_|_|\\__|_|___/\\___\\__,_|_|\\___|     |_/	\n"
			<< "\n"
			<< "                             M4_Multiscale\n"
			<< "                     Copyright (2017) by Wentao Xu\n"
			<< "\n"
			<< " M4 stands for Multiscale Mesh-based Meshfree Method\n"
			<< " M4 is an universal solver for particle systems:\n"
			<< "\n"
			<< " Any questions, please contact:\n"
			<< " Wentao Xu   wx2151@columbia.edu\n"
			<< "\n"
			<< "************************************************************************\n"
			<< "\n";
	}



	void MetaSolver::printInfo() const {
		std::cout
			<< "***************************  FILE LOCATION  ****************************\n"
			<< "Target folder: " << folder_
			<< "\t File name: " << filename_ << std::endl
			<< "\n";
		std::cout
			<< "************************  DOMAIN CONFIGURATION  ************************\n"
			<< domain_config_ << std::endl << "\n"
			<< "*************************  ENSEMBLE INFOMATION  ************************\n"
			<< "** Slave **\n"
			<< "Total particles: " << slave_solv_->np_ << std::endl
			<< "Inter-particle distance: " << slave_solv_->dp_ << std::endl
			<< "Pos min: " << (slave_solv_->ensemble_->getPosMin()).transpose().format(CleanFmt) << std::endl
			<< "Pos max: " << (slave_solv_->ensemble_->getPosMax()).transpose().format(CleanFmt) << std::endl
			<< "\n** Master **\n"
			<< "Total particles: " << master_solv_->np_ << std::endl
			<< "Inter-particle distance: " << master_solv_->dp_ << std::endl
			<< "Pos min: " << (master_solv_->ensemble_->getPosMin()).transpose().format(CleanFmt) << std::endl
			<< "Pos max: " << (master_solv_->ensemble_->getPosMax()).transpose().format(CleanFmt) << std::endl
			<< "\n";
		std::cout
			<< "************************  COMPUTATION DETAILS  *************************\n";
		std::cout
			<< "** Slave **\n";
		if (slave_solv_->l_compute_ != nullptr) {
			std::cout << "Lagrangian Solver: " << slave_solv_->l_compute_->getClassName() << std::endl;
			std::cout << "Lagrangian compute params: " << slave_solv_->l_compute_->format();
		}
		if (slave_solv_->e_compute_ != nullptr) {
			std::cout << "Eulerian solver: " << slave_solv_->e_compute_->getClassName() << std::endl;
			std::cout << "Lagrangian compute params: " << slave_solv_->e_compute_->format();
		}
		std::cout << std::endl;
		std::cout
			<< "\n** Master **\n";
		if (master_solv_->l_compute_ != nullptr) {
			std::cout << "Lagrangian Solver: " << master_solv_->l_compute_->getClassName() << std::endl;
			std::cout << "Lagrangian compute params: " << master_solv_->l_compute_->format();
		}
		if (master_solv_->e_compute_ != nullptr) {
			std::cout << "Eulerian solver: " << master_solv_->e_compute_->getClassName() << std::endl;
			std::cout << "Lagrangian compute params: " << master_solv_->e_compute_->format();
		}
		std::cout << std::endl;
	}

	void MetaSolver::printMemoryInfo() const {
		double s_nb = (slave_solv_->nbh_ != nullptr) ? slave_solv_->nbh_->getMemoryInMB() : 0.0;
		double s_pb = (slave_solv_->pbh_ != nullptr) ? slave_solv_->pbh_->getMemoryInMB() : 0.0;
		double m_nb = (master_solv_->nbh_ != nullptr) ? master_solv_->nbh_->getMemoryInMB() : 0.0;
		double m_pb = (master_solv_->pbh_ != nullptr) ? master_solv_->pbh_->getMemoryInMB() : 0.0;
		std::cout
			<< "***************************  MEMORY SUMMARY  ***************************\n"
			<< "** Slave **\n"
			<< "Allocated memory for particle: " << slave_solv_->ensemble_->getMemoryInMB() << " MB" << std::endl
			<< "Allocated memory for sorting:  " << slave_solv_->sorted_->getMemoryInMB() << " MB" << std::endl
			<< "Allocated memory for bonds: " << s_nb + s_pb << " MB" << std::endl
			<< "Total allocated memory: " << slave_solv_->ensemble_->getMemoryInMB() + slave_solv_->sorted_->getMemoryInMB() 
			+ s_nb + s_pb << " MB" << std::endl
			<< "\n";
		std::cout 
			<< "** Master **\n"
			<< "Allocated memory for particle: " << master_solv_->ensemble_->getMemoryInMB() << " MB" << std::endl
			<< "Allocated memory for sorting:  " << master_solv_->sorted_->getMemoryInMB() << " MB" << std::endl
			<< "Allocated memory for bonds: " << m_nb + m_pb << " MB" << std::endl
			<< "Total allocated memory: " << master_solv_->ensemble_->getMemoryInMB() + master_solv_->sorted_->getMemoryInMB()
			+ m_nb + m_pb << " MB" << std::endl
			<< "\n";
		std::cout
			<< "Total allocated memory for meta solver: " << slave_solv_->ensemble_->getMemoryInMB() +
			slave_solv_->sorted_->getMemoryInMB() + s_nb + s_pb + master_solv_->ensemble_->getMemoryInMB() +
			master_solv_->sorted_->getMemoryInMB() + m_nb + m_pb << " MB" << std::endl;
	}

	void MetaSolver::configMeshfreeSolver(std::string slav_lc, std::string slav_ec, std::vector<std::string> slav_attrs,
		std::string master_lc, std::string master_ec, std::vector<std::string> master_attrs) {
		std::cout
			<< "************************  CONFIGURING SOLVER  **************************\n"
			<< "configuring slave ..\n";
		slave_solv_->configSolver(slav_lc, slav_ec, slav_attrs);
		std::cout << "slave configuration done!\n";
		std::cout << "configuring master ..\n";
		master_solv_->configSolver(master_lc, master_ec, master_attrs);
		std::cout << "master configuration done!\n";
	}

	void MetaSolver::configContact() {
		if (slave_solv_ == nullptr || master_solv_ == nullptr)
			throwException("configContact", "uninitialized solver!");
		slave_sort_ = slave_solv_->sorted_;
		master_sort_ = master_solv_->sorted_;
		if (slave_sort_ == nullptr || master_sort_ == nullptr)
			throwException("configContact", "ensemble unsorted");
		if (master_solv_->l_compute_->getClassName() == "RigidBody")
		cm_ = std::make_shared<ContactManager>(master_sort_, slave_sort_, (dp_s_ + dp_m_) / 1.8, dp_s_, dt_);
		cm_->setForceParams(dp_s_ * dp_s_ * dp_s_, 10 * pow(10, 11)); // needs modification
	}

	void MetaSolver::configSRB() {
		use_srb_ = true;
		if (slave_solv_ == nullptr || master_solv_ == nullptr)
			throwException("configSRB", "uninitialized solver!");
		slave_sort_ = slave_solv_->sorted_;
		//master_sort_ = master_solv_->sorted_;
		if (slave_sort_ == nullptr)
			throwException("configSRB", "ensemble unsorted");
		if (master_solv_->l_compute_->getClassName() != "RigidBody")
			throwException("configSRB", "master not rigid body!");
		std::shared_ptr<RigidBody> rgb = std::static_pointer_cast<RigidBody>(master_solv_->l_compute_);
		double radius = master_obj_.dims[0];
		srb_cm_ = std::make_shared<ContactMngerSRB>(rgb, slave_sort_, radius, dp_s_, 5.0 * pow(10, 5));
	}

	void MetaSolver::createScene(ObjInfo master, ObjInfo slave, DomainConfig domain_config) {
		master_obj_ = master;
		slave_obj_ = slave;
		slave_solv_->createEnsemble(slave, domain_config);
		master_solv_->createEnsemble(master, domain_config);
		setDomainConfig(domain_config);
		dp_s_ = slave.dp;
		dp_m_ = master.dp;
	}

	void MetaSolver::configRun(double dt, double total_time, int sv_step) {
		sv_step_ = sv_step;
		dt_max_ = dt_ = dt;
		T_ = total_time;
		t_intv_ = sv_step_ * dt_;
		t_cur_ = 0.0;
		timestep_ = 0;
		sim_time_ = constant_timestep_ ? (int(T_ / dt_) + 1) * dt_ : (int(T_ / t_intv_) + 1) * t_intv_;
		slave_solv_->configRun(dt, total_time, sv_step);
		master_solv_->configRun(dt, total_time, sv_step);
	}

	void MetaSolver::addConfinedRegions(std::vector<std::shared_ptr<Shape>> regions) {
		for (int i = 0; i < regions.size(); i++)
			slave_solv_->configPosConstraint(regions[i]);
	}

	void MetaSolver::addConfinedRegionsComplement(std::vector<std::shared_ptr<Shape>> regions) {
		for (int i = 0; i < regions.size(); i++)
			slave_solv_->configPosConstraintComplement(regions[i]);
	}

	void MetaSolver::run() {
		printInfo();
		printMemoryInfo();
		std::cout
			<< "========================================================================\n"
			<< "-------------------------  M4 SOLVER RUNNING  --------------------------\n"
			<< "........................................................................\n";
		std::cout
			<< "Current time: " << theWatch().getDayTime()
			<< "    Date: " << theWatch().getDate() << "\n";
		sv_count_ = 0;
		if (constant_timestep_) {
			for (; timestep_ <= int(T_ / dt_); timestep_++) {
				if (timestep_ % sv_step_ == 0) {
					saveVtk();
					std::cout << " SaveVtk: " << sv_count_ << "\tTimestep: " << timestep_
						<< "\t" << std::endl;
				}
				theWatch().startTimer("run");
				verletUpdate();
				if (timestep_ % 10 == 0)
					std::cout << "(Constant) Timestep:\t" << timestep_ << "\t\t" << theWatch().stopAndFormat("run") << std::endl;
			}
			saveVtk();
			std::cout << " SaveVtk " << sv_count_ << " finished! " << std::endl;

		}
		else {
			int cnt = int(T_ / t_intv_);
			for (; timestep_ <= cnt; timestep_++) {
				saveVtk();
				std::cout << " SaveVtk: " << sv_count_ << "\tTimestep: " << timestep_
					<< "\t" << std::endl;
				while (t_cur_ < t_intv_) {
					theWatch().startTimer("run");
					verletUpdate();
					if (timestep_ % 10 == 0)
						std::cout << "(Variable) Timestep:\t" << timestep_  << " cur/intv: "
						<< t_cur_/t_intv_ 
						<< "\t\t" << theWatch().stopAndFormat("run") << std::endl;
					t_cur_ += cm_->getDt();
				}
				t_cur_ = 0.0;
			}
			saveVtk();
			std::cout << " SaveVtk: " << sv_count_ << "\tTimestep: " << timestep_
				<< "\t" << std::endl;
		}
	}

	void MetaSolver::runWithNoSlip() {
		printInfo();
		printMemoryInfo();
		std::cout
			<< "========================================================================\n"
			<< "-------------------------  M4 SOLVER RUNNING  --------------------------\n"
			<< "........................................................................\n";
		std::cout
			<< "Current time: " << theWatch().getDayTime()
			<< "    Date: " << theWatch().getDate() << "\n";
		sv_count_ = 0;
		if (constant_timestep_) {
			for (; timestep_ <= int(T_ / dt_); timestep_++) {
				if (timestep_ % sv_step_ == 0) {
					saveVtk();
					std::cout << " SaveVtk: " << sv_count_ << "\tTimestep: " << timestep_
						<< "\t" << std::endl;
				}
				theWatch().startTimer("run");
				verletUpdateWithNoSlip();
				if (timestep_ % 10 == 0)
					std::cout << "(Constant) Timestep:\t" << timestep_ << "\t\t" << theWatch().stopAndFormat("run") << std::endl;
			}
			saveVtk();
			std::cout << " SaveVtk " << sv_count_ << " finished! " << std::endl;

		}
		else {
			int cnt = int(T_ / t_intv_);
			for (; timestep_ <= cnt; timestep_++) {
				saveVtk();
				std::cout << " SaveVtk: " << sv_count_ << "\tTimestep: " << timestep_
					<< "\t" << std::endl;
				while (t_cur_ < t_intv_) {
					theWatch().startTimer("run");
					verletUpdateWithNoSlip();
					if (timestep_ % 10 == 0)
						std::cout << "(Variable) Timestep:\t" << timestep_ << " cur/intv: "
						<< t_cur_ / t_intv_
						<< "\t\t" << theWatch().stopAndFormat("run") << std::endl;
					t_cur_ += cm_->getDt();
				}
				t_cur_ = 0.0;
			}
			saveVtk();
			std::cout << " SaveVtk: " << sv_count_ << "\tTimestep: " << timestep_
				<< "\t" << std::endl;
		}
	}

	void MetaSolver::saveVtk() {
		saveVtkSlave();
		saveVtkMaster();
		saveVtkAll();
		sv_count_++;
	}

	void MetaSolver::saveVtkAll() {
		std::string inner_folder = folder_ + "\\" + filename_ + "_All";
		std::string vtkname = filename_ + "_All_" + std::to_string(sv_count_) + ".vtk";
		std::string damage_front_txt = filename_ + "_damage_front_" + std::to_string(sv_count_) + ".txt";
		std::ofstream damage_file(inner_folder + "\\" + damage_front_txt, std::ofstream::binary);
		for (int i = 0; i < 10; i++) {
			damage_file << slave_solv_->damage_front_[i] << " ";
		}
		damage_file.close();
		int s_np = slave_solv_->np_;
		int m_np = master_solv_->np_;
		long np = s_np + m_np;
		std::ofstream file(inner_folder + "\\" + vtkname, std::ios::binary);
		file << "# vtk DataFile Version 3.0\n"
			<< "vtk output\n"
			<< "ASCII\n"
			<< "DATASET POLYDATA\n"
			<< "POINTS " << np << " double\n";
		for (int i = 0; i < s_np; i++)
			file << slave_solv_->pos_[i].format(PureFmt) << "\n";
		for (int i = 0; i < m_np; i++)
			file << master_solv_->pos_[i].format(PureFmt) << "\n";
		file << "VERTICES " << np << " " << np * 2 << "\n";
		for (int i = 0; i < np; i++)
			file << "1 " << i << "\n";
		file << "POINT_DATA " << np << "\n";
		file << "SCALARS Object_ID int\n";
		file << "LOOKUP_TABLE default\n";
		for (int i = 0; i < s_np; i++)
			file << "0\n";
		for (int i = 0; i < m_np; i++)
			file << "1\n";

		file << "VECTORS Vel double\n";
		for (int i = 0; i < s_np; i++)
			file << slave_solv_->vel_[i].format(PureFmt) << "\n";
		for (int i = 0; i < m_np; i++)
			file << master_solv_->vel_[i].format(PureFmt) << "\n";

		file << "VECTORS Acc double\n";
		for (int i = 0; i < s_np; i++)
			file << slave_solv_->stack_[i].format(PureFmt) << "\n";
		for (int i = 0; i < m_np; i++)
			file << master_solv_->stack_[i].format(PureFmt) << "\n";
		file.close();
	}

	void MetaSolver::saveVtkSlave() {
		slave_solv_->saveVtk();
		slave_solv_->updateDamageFront();
	}

	void MetaSolver::saveVtkMaster() {
		master_solv_->saveVtk();
	}

	void MetaSolver::verletUpdate() {
		slave_sort_->makeSortPartial();
		master_sort_->makeSortPartial();
		//std::cout << "compute contact..\n";
		if (!use_srb_) {
			cm_->updateContactZone();
			cm_->computeContact();
			//std::cout << "compute contact done!\n";
			if (!constant_timestep_) {
				dt_ = cm_->getDt();
				slave_solv_->setDt(dt_);
				master_solv_->setDt(dt_);
			}
		}
		else {
			srb_cm_->computeForces();
		}
		//std::cout << "updating master..\n";
		
		//std::cout << "master updated!\n";
		//std::cout << "updating slave..\n";
		if (!use_srb_)
			slave_solv_->verletUpdate();
		else {
			slave_solv_->verletUpdateWithNoSlip(srb_cm_->getCenterPos(), srb_cm_->getCenterVel(), 
				srb_cm_->getCenterAcc(), srb_cm_->getEpsilon() - srb_cm_->getDp());
		}
		master_solv_->verletUpdate();
		//std::cout << "slave updated!\n"; 
		
	}

	void MetaSolver::verletUpdateWithNoSlip() {
		slave_sort_->makeSortPartial();
		master_sort_->makeSortPartial();
		slave_solv_->verletUpdateStepOne();
		master_solv_->verletUpdateStepOne();
		slave_solv_->verletComputeForces();
		master_solv_->verletComputeForces();
		srb_cm_->computeForces();
		slave_solv_->verletUpdateStepTwo();
		master_solv_->verletUpdateStepTwo();
		//-------------------
		//slave_solv_->verletUpdateConfineVel(srb_cm_->getDp() * 0.04 / dt_);
		double master_mass = srb_cm_->getMasterTotalMass();
		double eps = srb_cm_->getRadius()- srb_cm_->getDp();
		srb_cm_->update();
		Vec3d master_center = srb_cm_->getCenterPos();
		Vec3d vel = srb_cm_->getCenterVel();
		Vec3d acc = srb_cm_->getCenterAcc();
		Vec3d momentum = vel * master_mass;
		Vec3d impulse = acc * master_mass;
		slave_solv_->enforceNoSlip(momentum, impulse, master_mass, master_center, eps);
		slave_solv_->verletUpdateStepThree();
		//slave_solv_->verletUpdateConfineAcc(srb_cm_->getDp() * 0.9/ (dt_ * dt_));
		master_solv_->verletUpdateStepThree();
	}

	void MetaSolver::configIO(std::string folder, std::string filename, 
		std::string slave, std::string master) {
		slave_solv_->configIO(folder, filename + "_" + slave);
		master_solv_->configIO(folder, filename + "_" + master);
		folder_ = folder;
		filename_ = filename;
		std::string inner_folder = folder_ + "\\" + filename_ + "_All";
		CreateDirectory(std::wstring(inner_folder.begin(), inner_folder.end()).c_str(), NULL);
		slave_name_ = slave;
		master_name_ = master;		
	}
}