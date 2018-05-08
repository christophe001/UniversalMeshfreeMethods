/*! \file m_solver.cpp */

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

#include "m_solver.h"
#include "u_timer.h"
#include <windows.h>
#include <fstream>
#include <cmath>

constexpr auto theWatch = msl::Timer::instance;

namespace msl {
	ComputeSolver::ComputeSolver()
		: timestep_(0), max_step_(0), sv_step_(50),
		sv_count_(0), compute_(0)
	{
		class_name_ = "ComputeSolver";
	}

	void ComputeSolver::printLabel() {
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

	void ComputeSolver::printInfo() {
		std::cout
			<< "***************************  FILE LOCATION  ****************************\n"
			<< "Target folder: " << folder_
			<< "\t File name: " << filename_ << std::endl
			<< "\n";
		std::cout
			<< "************************  DOMAIN CONFIGURATION  ************************\n"
			<< *(sorted_->getDomainConfig()) << std::endl
			<< "*************************  ENSEMBLE INFOMATION  ************************\n"
			<< "Total particles: " << np_ << std::endl
			<< "Inter-particle distance: " << dp_ << std::endl
			<< "\n";
	}

	void ComputeSolver::printMemoryInfo() {
		std::cout
			<< "***************************  MEMORY SUMMARY  ***************************\n"
			<< "Allocated memory for particle: " << ensemble_->getMemoryInMB() << std::endl
			<< "Allocated memory for sorting:  " << sorted_->getMemoryInMB() << std::endl
			<< "Allocated memory for bonds: " << nbh_->getMemoryInMB() + pbh_->getMemoryInMB() << std::endl
			<< "Allocated memory for computing: " << compute_->getMemoryInMB() << std::endl
			<< "Total memory: " << ensemble_->getMemoryInMB() + sorted_->getMemoryInMB() + nbh_->getMemoryInMB()
			+ pbh_->getMemoryInMB() + compute_->getMemoryInMB() << std::endl
			<< "\n";
	}

	void ComputeSolver::create2DCasePeridynamics(Peridm2DParams pp, bool ec) {
		std::cout
			<< "========================================================================\n"
			<< "-------------------------  Creating Ensemble  --------------------------\n"
			<< "........................................................................\n";
		theWatch().startTimer("create");
		ensemble_ = std::make_shared<Ensemble>();
		ensemble_->addVectorAttribute(0);
		ensemble_->addScalarAttribute("Initial_position");
		ensemble_->addScalarAttribute("Damage_ratio");
		horizon_ = pp.horizon;
		dp_ = horizon_ / pp.m;
		double length_ = pp.length;
		double width_ = pp.width;
		int nl = width_ / (2.0 * dp_);
		dp_ = width_ / (2.0 * nl);
		nl = nl * 2;
		int nw = length_ / dp_;
		length_ = nw * dp_;
		np_ = nw * nl;
		ensemble_->setNp(np_);
		Vec3d offset{ (length_ - dp_) / 2.0, 0.0, (width_ - dp_) / 2.0 };
		DomainConfig cfg(-2.0 * offset, 2.0 * offset, horizon_, PeriType::kPeriO);
		ensemble_->setDomainConfig(cfg);
		pos_ = ensemble_->getPos();
		vel_ = ensemble_->getVel();
		acc_ = ensemble_->getAcc();
		scalar_attrs_ = ensemble_->getScalarAttrPtrs();
		vector_attrs_ = ensemble_->getVectorAttrPtrs();
		Vec3d* ipos = vector_attrs_[0]->getAttr();
#pragma omp parallel for schedule(static)
		for (int i = 0; i < np_; i++) {
			pos_[i] = Vec3d{ (i % nw) * dp_, 0, (i / nw) * dp_ } -offset;
			ipos[i] = pos_[i];
		}
		std::cout
			<< "Setting particles: " << np_ << " " << theWatch().stopAndFormat("create") << std::endl;
		theWatch().startTimer("sort");
		sorted_ = std::make_shared<SortEnsemble>(ensemble_);
		sorted_->makeSortFull();
		std::cout
			<< "Sort particles:    " << theWatch().stopAndFormat("sort") << std::endl;
		theWatch().startTimer("bond");
		std::cout
			<< "Total cells: " << sorted_->getTotalCells() << std::endl;
		compute_bw_ = std::make_shared<Peridm2D>(sorted_, pp, ec);
		compute_bw_->initialize();
		std::cout
			<< "Computing bond:    " << theWatch().stopAndFormat("bond") << std::endl;
		nbh_ = compute_bw_->getNeighborhoodData();
		pbh_ = compute_bw_->getPeriNeighborData();

	}

	double ComputeSolver::getModulus() {
		double maxacc = getMaxAcc();
		Vec3d z1 = { 0,0,1 };
		Vec3d z2 = { 0,0,-1 };
		double h1 = getMaxPos(z1).dot(z1) - getMaxPos(z2).dot(z1);
		std::cout << "initial width: " << h1 << std::endl;
		int t = 0;
		while (maxacc > pow(10, -4)) {
			verletUpdate();
			maxacc = getMaxAcc();
			double ht = getMaxPos(z1).dot(z1) - getMaxPos(z2).dot(z1);
			if (t % 20 == 0) {
				std::cout << "width: " << ht << std::endl;
				std::cout << "strain " << (h1 - ht) / h1 << std::endl;
			}
			t++;
		}
		double h2 = getMaxPos(z1).dot(z1) - getMaxPos(z2).dot(z1);
		std::cout << "final width: " << h2 << std::endl;
		if (h1 == h2)
			return 0;
		else return h2 / h1 - 1.0;
	}

	double ComputeSolver::getMaxAcc() const {
		double ans = 0;
		for (int i = 0; i < np_; i++) {
			ans = max(ans, acc_[i].norm());
		}
		return ans;
	}

	double ComputeSolver::getMaxVel() const {
		double ans = 0;
		for (int i = 0; i < np_; i++) {
			ans = max(ans, vel_[i].norm());
		}
		return ans;
	}

	Vec3d ComputeSolver::getMaxPos(Vec3i dir) const {
		Vec3d ans = pos_[0];
		double t = ans.dot(dir);
		double temp;
		for (int i = 1; i < np_; i++) {
			temp = pos_[i].dot(dir);
			if (t < temp) {
				t = temp;
				ans = pos_[i];
			}
		}
		return ans;
	}

	void ComputeSolver::run() {
		if (!compute_ && !compute_bw_) return;
		printLabel();
		printInfo();
		std::cout
			<< "========================================================================\n"
			<< "-------------------------  M4 SOLVER RUNNING  --------------------------\n"
			<< "........................................................................\n";
		std::cout
			<< "Current time: " << theWatch().getCurrentTime()
			<< " Date: " << theWatch().getDayTime() << "\n";

		for (; timestep_ <= max_step_; timestep_++) {
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
		if (timestep_ % sv_step_ != 0) {
			saveVtk();
			std::cout << " SaveVtk " << sv_count_ << std::endl;
		}
	}

	void ComputeSolver::verletUpdate() {
		if (compute_bw_ != 0)
			compute_bw_->verletUpdate(dt_);
		else compute_->verletUpdate(dt_);
	}

	void ComputeSolver::configParams(std::string folder, std::string filename, double dt, double total_time) {
		if (!ensemble_)
			throwException("configParams", "cannot find underlying particle ensemble");
		folder_ = folder;
		filename_ = filename;
		dt_ = dt;
		T_ = int(total_time / dt) * dt;
		timestep_ = 0;
		max_step_ = int(total_time / dt);
		std::string inner_folder = folder_ + "\\" + filename_;
		CreateDirectory(std::wstring(folder_.begin(), folder_.end()).c_str(), NULL);
		CreateDirectory(std::wstring(inner_folder.begin(), inner_folder.end()).c_str(), NULL);
	}

	void ComputeSolver::saveVtk() {
		std::string inner_folder = folder_ + "\\" + filename_;
		std::string vtkname = filename_ + std::to_string(sv_count_) + ".vtk";
		double sim_time = timestep_ * dt_;

		std::ofstream file(inner_folder + "\\" + vtkname, std::ios::binary);
		file << "# vtk DataFile Version 3.0 \n"
			<< "vtk output \n"
			<< "ASCII \n"
			<< "DATASET POLYDATA \n"
			<< "POINTS " << np_ << " double\n";
		for (int i = 0; i < np_; i++) {
			file << pos_[i].format(PureFmt) << "\n";
		}
		file << "VERTICES " << np_ << " " << np_ * 2 << "\n";
		for (int i = 0; i < np_; i++)
			file << "1 " << i << " \n";
		file << "POINT_DATA " << np_ << " \n";
		for (auto&& sa : scalar_attrs_) {
			file << "SCALARS " << sa->getInfo() << " double\n";
			file << "LOOKUP_TABLE default\n";
			for (int i = 0; i < np_; i++)
				file << sa->getAttr()[i] << "\n";
		}

		file << "VECTORS vel double\n";
		for (int i = 0; i < np_; i++) {
			file << vel_[i].format(PureFmt) << "\n";
		}

		file << "VECTORS acc double\n";
		for (int i = 0; i < np_; i++) {
			file << acc_[i].format(PureFmt) << " \n";
		}

		for (auto&& va : vector_attrs_) {
			file << "VECTORS " << va->getInfo() << " double\n";
			for (int i = 0; i < np_; i++)
				file << va->getAttr()[i].format(PureFmt) << " \n";
		}
		file.close();
		sv_count_++;
	}
}