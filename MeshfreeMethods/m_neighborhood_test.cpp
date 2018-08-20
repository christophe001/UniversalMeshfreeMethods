/*! \file m_neighborhood_test.cpp */

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

#include "m_neighborhood_test.h"
#include <windows.h>
#include <fstream>
#include <iostream>
#include "u_timer.h"
#include <random>

namespace msl {
	NeighborhoodTest::NeighborhoodTest(bool sim2d, Vec3d dims, double dp) {
		sim2d_ = sim2d;
		dims_ = dims;
		if (sim2d) dims_[1] = 0;
		dp_ = dp;
		class_name_ = "NeighborhoodTest";
	}

	void NeighborhoodTest::run() {
		constexpr auto theWatch = msl::Timer::instance;
		std::cout << "creating ensemble...\n";
		theWatch().startTimer("create");
		if (sim2d_) {
			creator_ = std::make_shared<EnsembleCreator>("rectangle2d", Vec3d::Zero());
		}
		else {
			creator_ = std::make_shared<EnsembleCreator>("rectangle", Vec3d::Zero());
		}
		std::vector<std::string> sattr;
		sattr.push_back("within");
		creator_->addScalarAttributes(sattr);
		creator_->setDims(dims_);
		creator_->setDpActual(dp_);
		creator_->create();
		creator_->setScalarAttributeConstant("within", 0.1);
		ensemble_ = creator_->getEnsemble();
		np_ = ensemble_->getNp();
		pos_ = ensemble_->getPos();
		vel_ = ensemble_->getVel();
		acc_ = ensemble_->getAcc();
		std::cout << "create done!    total time: " << theWatch().stopAndFormat("create") << std::endl;
		std::cout << "total particles: " << ensemble_->getNp() << std::endl;
		std::cout << "scalar attributes: ";
		for (auto sa : ensemble_->getScalarAttrPtrs()) {
			std::cout << sa->getInfo() << " ";
		}
		std::cout << std::endl;
		//! sort particles
		theWatch().startTimer("sort");
		std::cout << "\n sorting particles...\n";
		cpn_ = std::make_shared<ComputeNeighbor>(ensemble_);
		Vec3d pmin{ -10, -10, -10 };
		Vec3d pmax{ 10, 10, 10 };
		if (creator_->getActualShape()->shape2D()) {
			pmin[1] = 0;
			pmax[1] = 0;
		}

		DomainConfig domain_config(pmin, pmax, 1.0);
		cpn_->setDomainConfig(domain_config);
		cpn_->makeSortFull();
		nbh_ = cpn_->getNeighborhoodData();
		pbh_ = cpn_->getPeriNeighborData();
		std::cout << "sort particles done!    total time: " << theWatch().stopAndFormat("sort") << std::endl;
		//! compute bond
		theWatch().startTimer("compute");
		std::cout << "\n computing bond...\n";
		std::vector<std::string> attr;
		cpn_->initialize(dp_, horizon_, attr);
		cpn_->compute();
		std::cout << "compute bond done!    total time: " << theWatch().stopAndFormat("compute") << std::endl;
		std::cout << "total bonds: " << cpn_->getNeighborhoodData()->getTotalBonds() << std::endl;
		//! add crack
		theWatch().startTimer("crack");
		std::cout << "\n setting crack...\n";
		cpn_->addSplit(Vec3d::UnitZ(), 0, Vec3d::UnitX(), 0);
		cpn_->updateDamage();
		std::cout << "set crack done!    total time: " << theWatch().stopAndFormat("crack") << std::endl;
		//! within
		std::uniform_int_distribution<> dis(0, np_ - 1);
		std::default_random_engine re;
		re.seed(std::chrono::system_clock::now().time_since_epoch().count());
		double* within = ensemble_->getScalarAttrPtr("within")->getAttr();
		for (int i = 0; i < 6; i++) {
			int id = dis(re);
			std::cout << "particle id: " << id << std::endl;
			within[id] = 1.0;
			for (long it = nbh_->begin(id); it != nbh_->end(id); it++) {
				within[nbh_->getNeighborhoodList()[it]] = 1;
			}
		}

	}

	void NeighborhoodTest::runPeri() {
		constexpr auto theWatch = msl::Timer::instance;
		std::cout << "creating ensemble...\n";
		theWatch().startTimer("create");
		if (sim2d_) {
			creator_ = std::make_shared<EnsembleCreator>("rectangle2d", Vec3d::Zero());
		}
		else {
			creator_ = std::make_shared<EnsembleCreator>("rectangle", Vec3d::Zero());
		}
		std::vector<std::string> sattr;
		sattr.push_back("within");
		creator_->addScalarAttributes(sattr);
		creator_->setDims(dims_);
		creator_->setDpActual(dp_);
		creator_->create();
		creator_->setScalarAttributeConstant("within", 0.1);
		ensemble_ = creator_->getEnsemble();
		np_ = ensemble_->getNp();
		pos_ = ensemble_->getPos();
		vel_ = ensemble_->getVel();
		acc_ = ensemble_->getAcc();
		std::cout << "create done!    total time: " << theWatch().stopAndFormat("create") << std::endl;
		std::cout << "total particles: " << ensemble_->getNp() << std::endl;
		std::cout << "scalar attributes: ";
		for (auto sa : ensemble_->getScalarAttrPtrs()) {
			std::cout << sa->getInfo() << " ";
		}
		std::cout << std::endl;
		//! sort particles
		theWatch().startTimer("sort");
		std::cout << "\n sorting particles...\n";
		cpn_ = std::make_shared<ComputeNeighbor>(ensemble_);
		Vec3d pmin = -0.5 * dims_;
		Vec3d pmax = 0.5 * dims_;
		DomainConfig domain_config(pmin, pmax, 1.0, PeriType::kPeriXYZ);
		cpn_->setDomainConfig(domain_config);
		cpn_->makeSortFull();
		nbh_ = cpn_->getNeighborhoodData();
		pbh_ = cpn_->getPeriNeighborData();
		std::cout << "sort particles done!    total time: " << theWatch().stopAndFormat("sort") << std::endl;
		//! compute bond
		theWatch().startTimer("compute");
		std::cout << "\n computing bond...\n";
		std::vector<std::string> attr;
		cpn_->initialize(dp_, horizon_, attr);
		cpn_->compute();
		std::cout << "compute bond done!    total time: " << theWatch().stopAndFormat("compute") << std::endl;
		std::cout << "total bonds: " << cpn_->getNeighborhoodData()->getTotalBonds() << std::endl;
		//! add crack
		theWatch().startTimer("crack");
		std::cout << "\n setting crack...\n";
		cpn_->addSplit(Vec3d::UnitZ(), 0, Vec3d::UnitX(), 0);
		cpn_->updateDamage();
		std::cout << "set crack done!    total time: " << theWatch().stopAndFormat("crack") << std::endl;
		//! within
		std::uniform_int_distribution<> dis(0, np_ - 1);
		std::default_random_engine re;
		re.seed(std::chrono::system_clock::now().time_since_epoch().count());
		double* within = ensemble_->getScalarAttrPtr("within")->getAttr();
		for (int i = 0; i < 6; i++) {
			int id = (i == 0) ? 0 : dis(re);
			std::cout << "particle id: " << id << std::endl;
			within[id] = 1;
			for (long it = nbh_->begin(id); it != nbh_->end(id); it++) {
				within[nbh_->getNeighborhoodList()[it]] = 1;
			}
			if (pbh_->hasPeriBond(id)) {
				for (long it = pbh_->begin(id); it != pbh_->end(id); it++) {
					within[pbh_->getPeriBondList()[it].m_j] = 0.5;
				}
			}
		}
	}

	void NeighborhoodTest::setSaveParams(std::string folder, std::string file) {
		folder_ = folder;
		file_ = file;
		std::string inner_folder = folder_ + "\\" + file_;
		CreateDirectory(std::wstring(folder_.begin(), folder_.end()).c_str(), NULL);
		CreateDirectory(std::wstring(inner_folder.begin(), inner_folder.end()).c_str(), NULL);
	}

	void NeighborhoodTest::saveVtk() {
		std::string inner_folder = folder_ + "\\" + file_;
		std::string vtkname = file_ + ".vtk";
		std::ofstream file(inner_folder + "\\" + vtkname, std::ios::binary);
		file << "# vtk DataFile Version 3.0\n"
			<< "vtk output\n"
			<< "ASCII\n"
			<< "DATASET POLYDATA\n"
			<< "POINTS " << np_ << " double\n";
		for (int i = 0; i < np_; i++) {
			file << pos_[i].format(PureFmt) << "\n";
		}
		file << "VERTICES " << np_ << " " << np_ * 2 << "\n";
		for (int i = 0; i < np_; i++)
			file << "1 " << i << "\n";
		file << "POINT_DATA " << np_ << "\n";
		int* bc = nbh_->getBondCount();
		double* damage = nbh_->getParticleDamage();
		file << "SCALARS damage double\n";
		file << "LOOKUP_TABLE default\n";
		for (int i = 0; i < np_; i++)
			file << damage[i]/(double)bc[i] << "\n";

		double* within = ensemble_->getScalarAttrPtr("within")->getAttr();
		file << "SCALARS within double\n";
		file << "LOOKUP_TABLE default\n";
		for (int i = 0; i < np_; i++)
			file << within[i] << "\n";
		file << "VECTORS Vel double\n";
		for (int i = 0; i < np_; i++) {
			file << vel_[i].format(PureFmt) << "\n";
		}

		file << "VECTORS Acc double\n";
		for (int i = 0; i < np_; i++) {
			file << acc_[i].format(PureFmt) << "\n";
		}
		file.close();
	}
}