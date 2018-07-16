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
	MetaSolver::MetaSolver() {
		class_name_ = "MetaSolver";
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
			<< "\n** Master **\n"
			<< "Total particles: " << master_solv_->np_ << std::endl
			<< "Inter-particle distance: " << master_solv_->dp_ << std::endl
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
	}

	void MetaSolver::configRun(double dt, double total_time, int sv_step) {
		sv_step_ = sv_step;
		dt_max_ = dt_ = dt;
		T_ = total_time;
		t_intv_ = sv_step_ * dt_;
		t_cur_ = 0.0;
		timestep_ = 0;
		sim_time_ = constant_timestep_ ? (int(T_ / dt_) + 1) * dt_ : (int(T_ / t_intv_) + 1) * t_intv_;
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
					std::cout << "Timestep:\t" << timestep_ << "\t\t" << theWatch().stopAndFormat("run") << std::endl;
			}
			saveVtk();
			std::cout << " SaveVtk " << sv_count_ << std::endl;

		}
		else {
			int cnt = int(T_ / t_intv_);
			for (int i = 0; i <= cnt; i++) {
				saveVtk();
				std::cout << " SaveVtk: " << sv_count_ << "\tTimestep: " << timestep_
					<< "\t" << std::endl;
				while (t_cur_ < t_intv_) {
					theWatch().startTimer("run");
					verletUpdate();
					if (timestep_ % 10 == 0)
						std::cout << "Timestep:\t" << timestep_ << "\t\t" << theWatch().stopAndFormat("run") << std::endl;
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
	}

	void MetaSolver::saveVtkAll() {
		std::string inner_folder = folder_ + "\\" + filename_;
		std::string vtkname = filename_ + "_All_" + std::to_string(sv_count_) + ".vtk";
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
			file << slave_solv_->acc_[i].format(PureFmt) << "\n";
		for (int i = 0; i < m_np; i++)
			file << master_solv_->acc_[i].format(PureFmt) << "\n";
	}

	void MetaSolver::saveVtkSlave() {
		slave_solv_->saveVtk();
	}

	void MetaSolver::saveVtkMaster() {
		master_solv_->saveVtk();
	}

	void MetaSolver::verletUpdate() {
		cm_->updateContactZone();
		cm_->computeContact();
		if (!constant_timestep_) {
			dt_ = cm_->getDt();
			slave_solv_->setDt(dt_);
			master_solv_->setDt(dt_);
		}
		slave_solv_->verletUpdate();
		master_solv_->verletUpdate();
	}
	void MetaSolver::configIO(std::string folder, std::string filename, 
		std::string slave, std::string master) {
		slave_solv_->configIO(folder, filename + "_" + slave);
		master_solv_->configIO(folder, filename + "_" + master);
		folder_ = folder;
		filename_ = filename;
		std::string inner_folder = folder_ + "\\" + filename_;
		slave_name_ = slave;
		master_name_ = master;
		s_scalar_saves_ = slave_solv_->scalar_saves_;
		s_vector_saves_ = slave_solv_->vector_saves_;
		m_scalar_saves_ = master_solv_->scalar_saves_;
		m_vector_saves_ = master_solv_->vector_saves_;
		CreateDirectory(std::wstring(folder_.begin(), folder_.end()).c_str(), NULL);
		CreateDirectory(std::wstring(inner_folder.begin(), inner_folder.end()).c_str(), NULL);
	}
}