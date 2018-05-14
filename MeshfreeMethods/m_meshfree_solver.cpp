/*! \file m_meshfree_solver.cpp */

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

#include "m_meshfree_solver.h"
#include "u_timer.h"
#include <windows.h>
#include <fstream>
#include <cmath>

constexpr auto theWatch = msl::Timer::instance;

namespace msl {
	MeshfreeSolver::MeshfreeSolver() 
		: ensemble_(0), sorted_(0), nbh_(0), pbh_(0),
		l_compute_(0), e_compute_(0), creator_(0),
		id_(0), dict_(0), pos_(0), vel_(0), acc_(0),
		sv_count_(0), sv_step_(20)
	{
		relaxation_ = false;
		class_name_ = "MeshfreeSolver";
	}

	std::vector<MeshfreeSolver::AttributeInfo> MeshfreeSolver::getJH2AttributeConfig(double lambda, double mu) {
		auto attrs = std::vector<AttributeInfo>();
		typedef MeshfreeSolver::AttributeInfo MyAttr;
		attrs.push_back(MyAttr("tensor", "deformation", false, "identity"));
		attrs.push_back(MyAttr("tensor", "shape_tensor", false, "zero"));
		attrs.push_back(MyAttr("tensor", "deformation_dot", false, "zero"));
		attrs.push_back(MyAttr("tensor", "rotation", false, "zero"));
		attrs.push_back(MyAttr("tensor", "left_stretch", false, "zero"));
		attrs.push_back(MyAttr("tensor", "left_stretch", false, "zero"));
		attrs.push_back(MyAttr("tensor", "left_stretch_dot", false, "zero"));
		attrs.push_back(MyAttr("tensor", "d", false, "zero"));
		attrs.push_back(MyAttr("tensor", "tau", false, "zero"));
		attrs.push_back(MyAttr("vector", "initial_position", true, "pos"));
		attrs.push_back(MyAttr("scalar", "density", true, "constant", 3700));
		attrs.push_back(MyAttr("scalar", "lame_lambda", false, "constant", lambda));
		attrs.push_back(MyAttr("scalar", "lame_mu", false, "constant", mu));
		attrs.push_back(MyAttr("scalar", "pressure", true, "zero"));
		attrs.push_back(MyAttr("scalar", "damage", true, "zero"));
		attrs.push_back(MyAttr("scalar", "delta_p", false, "zero"));
		return attrs;
	}

	MeshfreeSolver::~MeshfreeSolver() {
		if (stack_ != 0) delete[] stack_;
	}

	void MeshfreeSolver::setParams(double horizon, double sv_step, bool relaxation, double rc) {
		horizon_ = horizon;
		sv_step_ = sv_step;
		relaxation_ = relaxation;
		rc_ = rc;
	}

	void MeshfreeSolver::setJH2Params(int i) {
		if (l_compute_ == 0)
			throwException("setJH2Params", "no lagrangian compute found");
		if (l_compute_->getClassName() != "StateBasedPDJH2")
			throwException("setJH2Params", "found " + l_compute_->getClassName() + " instead of JH2 compute");
		auto jh2_compute = std::static_pointer_cast<StateBasedPDJH2>(l_compute_);
		switch (i) {
		case 1:	l_compute_->configModelParams(jh2_param_1);
			break;
		case 2: l_compute_->configModelParams(jh2_param_2);
			break;
		case 3: l_compute_->configModelParams(jh2_param_3);
			break;
		default: throwException("setJH2Params", "unknown identifier");
		}
	}



	void MeshfreeSolver::printLabel() {
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

	void MeshfreeSolver::printInfo() {
		std::cout
			<< "***************************  FILE LOCATION  ****************************\n"
			<< "Target folder: " << folder_
			<< "\t File name: " << filename_ << std::endl
			<< "\n";
		std::cout
			<< "************************  DOMAIN CONFIGURATION  ************************\n"
			<< *(sorted_->getDomainConfig()) << std::endl << "\n"
			<< "*************************  ENSEMBLE INFOMATION  ************************\n"
			<< "Total particles: " << np_ << std::endl
			<< "Inter-particle distance: " << dp_ << std::endl
			<< "\n";
		std::cout
			<< "************************  COMPUTATION DETAILS  *************************\n";
		if (l_compute_ != 0) {
			std::cout << "Lagrangian Solver: " << l_compute_->getClassName() << std::endl;
			std::cout << "Lagrangian compute params: " << l_compute_->format();
		}
		if (e_compute_ != 0) {
			std::cout << "Eulerian solver: " << e_compute_->getClassName() << std::endl;
			std::cout << "Lagrangian compute params: " << e_compute_->format();
		}
		std::cout << std::endl;

	}

	void MeshfreeSolver::createEnsemble(std::string shape, DomainConfig domain_config, 
		std::vector<MeshfreeSolver::AttributeInfo>& infos, Vec3d dims, 
		Vec3d offset, Vec3d orientation, double theta, double dp) {
		printLabel();
		std::cout
			<< "*************************  CREATING ENSEMBLE  **************************\n"
			<< "creating..\n";
		creator_ = std::make_shared<EnsembleCreator>(shape, offset);
		std::vector<std::string> sas;
		std::vector<std::string> vas;
		std::vector<std::string> tas;
		for (auto info : infos) {
			if (info.m_attr_type == "scalar") {
				sas.push_back(info.m_name);
				if (info.m_save)
					scalar_saves_.push_back(info.m_name);
			}
			else if (info.m_attr_type == "vector") {
				vas.push_back(info.m_name);
				if (info.m_save)
					vector_saves_.push_back(info.m_name);
			}
			else if (info.m_attr_type == "tensor")
				tas.push_back(info.m_name);
			else throwException("createEnsemble", "Invalid attribute type!");
		}
		creator_->addScalarAttributes(sas);
		creator_->addVectorAttributes(vas);
		creator_->addTensorAttributes(tas);
		creator_->setDims(dims);
		creator_->setOrientation(orientation, theta);
		creator_->setDp(dp);
		creator_->create();
		for (auto info : infos) {
			if (info.m_attr_type == "scalar") {
				if (info.m_init == "zero")
					creator_->setScalarAttributeConstant(info.m_name, 0);
				else if (info.m_init == "one")
					creator_->setScalarAttributeConstant(info.m_name, 0);
				else if (info.m_init == "constant")
					creator_->setScalarAttributeConstant(info.m_name, info.m_val);
				else throwException("creatorEnsemble", "Invalid scalar initial condition!");
			}
			if (info.m_attr_type == "vector") {
				if (info.m_init == "zero")
					creator_->setVectorAttributeZero(info.m_name);
				else if (info.m_init == "pos") {
					Vec3d* ps = creator_->getEnsemble()->getPos();
					creator_->setVectorAttribute(info.m_name, ps);
				}
				else throwException("creatorEnsemble", "Invalid vector initial condition!");
			}
			if (info.m_attr_type == "tensor") {
				if (info.m_init == "zero")
					creator_->setTensorAttributeZero(info.m_name);
				else if (info.m_init == "identity")
					creator_->setTensorAttributeIdentity(info.m_name);
				else throwException("creatorEnsemble", "Invalid tensor initial condition!");
			}
		}
		ensemble_ = creator_->getEnsemble();
		pos_ = ensemble_->getPos();
		vel_ = ensemble_->getVel();
		acc_ = ensemble_->getAcc();
		np_ = ensemble_->getNp();
		id_ = ensemble_->getId();
		dict_ = ensemble_->getDict();
		dp_ = dp;
		try {
			stack_ = new Vec3d[np_];
		}
		catch (std::bad_alloc) {
			throwException("createEnsemble", "Error occured while allocating memory for stack");
		}
		scalar_attrs_ = ensemble_->getScalarAttrPtrs();
		vector_attrs_ = ensemble_->getVectorAttrPtrs();
		tensor_attrs_ = ensemble_->getTensorAttrPtrs();
		std::cout << "creating completed!\n";
		std::cout << "sorting ..\n";
		setDomainConfig(domain_config);
		sorted_ = std::make_shared<SortEnsemble>(ensemble_);
		sorted_->setDomainConfig(domain_config_);
		sorted_->makeSortFull(true);
		std::cout << "sorting completed!\n\n";

	}

	void MeshfreeSolver::configSolver(std::string lcompute, std::string ecompute,
		std::vector<std::string> attrs) {
		if (!lcompute.empty()) {
			std::cout
				<< "**********************  COMPUTING NEIGHBORHOOD  ************************\n"
				<< "computing ..\n";
			std::shared_ptr<ComputeNeighbor> cpn = std::make_shared<ComputeNeighbor>(sorted_);
			//std::cout << "dp: " << dp_ << " horizon: " << horizon_ << std::endl;
			cpn->initialize(dp_, horizon_, attrs);
			cpn->compute();
			nbh_ = cpn->getNeighborhoodData();
			pbh_ = cpn->getPeriNeighborData();
			std::cout << "neighborhood computation done!\n\n";
		}
		std::cout
			<< "***************************  CONFIG SOLVER  ****************************\n"
			<< "configuring ..\n";
		if (lcompute == "pd_jh2")
			l_compute_ = std::make_shared<StateBasedPDJH2>(sorted_, nbh_, pbh_);
		std::cout << "solver configuration succeed!\n\n";

	}

	void MeshfreeSolver::adaptiveRelaxation() {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++)
			acc_[i] -= rc_ * vel_[i];
	}

	void MeshfreeSolver::addPosForce(std::shared_ptr<Shape> shape, Vec3d force) {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++) {
			if (shape->isWithin(pos_[i])) {
				acc_[i] += force;
			}
		}
	}

	void MeshfreeSolver::addInitialPosForce(std::shared_ptr<Shape> shape, Vec3d force) {
		Vec3d * ipos = ensemble_->getVectorAttrPtr("initial_position")->getAttr();
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++) {
			if (shape->isWithin(ipos[i])) {
				acc_[dict_[i]] += force;
			}
		}
	}

	void MeshfreeSolver::configExternalPosForce(std::shared_ptr<Shape> shape, Vec3d force) {
		force_regions_.push_back(shape);
		forces_.push_back(force);
	}

	void MeshfreeSolver::configExternalInitPosForce(std::shared_ptr<Shape> shape, Vec3d force) {
		force_init_regions_.push_back(shape);
		forces_init_.push_back(force);
	}

	void MeshfreeSolver::initForces() {
		for (int i = 0; i < force_regions_.size(); i++) {
			addPosForce(force_regions_[i], forces_[i]);
		}
		for (int i = 0; i < force_init_regions_.size(); i++) {
			addInitialPosForce(force_init_regions_[i], forces_init_[i]);
		}
	}

	void MeshfreeSolver::printMemoryInfo() {
		std::cout
			<< "***************************  MEMORY SUMMARY  ***************************\n"
			<< "Allocated memory for particle: " << ensemble_->getMemoryInMB() << " MB" << std::endl
			<< "Allocated memory for sorting:  " << sorted_->getMemoryInMB() << " MB" << std::endl
			<< "Allocated memory for bonds: " << nbh_->getMemoryInMB() + pbh_->getMemoryInMB() << " MB" << std::endl
			<< "Total allocated memory: " << ensemble_->getMemoryInMB() + sorted_->getMemoryInMB() + nbh_->getMemoryInMB()
			+ pbh_->getMemoryInMB() << " MB" << std::endl
			<< "\n";
	}

	void MeshfreeSolver::verletUpdate() {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
			for (int i = 0; i < np_; i++) {
				stack_[i] = acc_[i];
				pos_[i] += vel_[i] * dt_ + 0.5 * dt_ * dt_ * acc_[i];
			}
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
			for (int i = 0; i < np_; i++)
				acc_[i] = Vec3d::Zero();
			if (l_compute_ != nullptr)
				l_compute_->computeForces();
			if (e_compute_ != nullptr)
				e_compute_->computeForces();
			if (relaxation_)
				adaptiveRelaxation();
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
			for (int i = 0; i < np_; i++)
				vel_[i] += (acc_[i] + stack_[i]) * 0.5 * dt_;
	}

	void MeshfreeSolver::saveVtk() {
		std::string inner_folder = folder_ + "\\" + filename_;
		std::string vtkname = filename_ + std::to_string(sv_count_) + ".vtk";
		double sim_time = timestep_ * dt_;

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
		for (auto&& s : scalar_saves_) {
			auto sa = ensemble_->getScalarAttrPtr(s);
			file << "SCALARS " << s << " double\n";
			file << "LOOKUP_TABLE default\n";
			for (int i = 0; i < np_; i++)
				file << sa->getAttr()[i] << "\n";		
		}

		file << "VECTORS Vel double\n";
		for (int i = 0; i < np_; i++) {
			file << vel_[i].format(PureFmt) << "\n";
		}

		file << "VECTORS Acc double\n";
		for (int i = 0; i < np_; i++) {
			file << acc_[i].format(PureFmt) << "\n";
		}

		for (auto&& s : vector_saves_) {
			auto va = ensemble_->getVectorAttrPtr(s);
			file << "VECTORS " << va->getInfo() << " double\n";
			for (int i = 0; i < np_; i++)
				file << va->getAttr()[i].format(PureFmt) << "\n";
		}
		file.close();
		sv_count_++;
	}

	void MeshfreeSolver::configRunParams(std::string folder, std::string filename, double dt, double total_time) {
		if (!ensemble_)
			throwException("configRunParams", "cannot find underlying particle ensemble");
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

	void MeshfreeSolver::run() {
		if (l_compute_ == 0 && e_compute_ == 0) return;
		//printLabel();
		printInfo();
		printMemoryInfo();
		std::cout
			<< "========================================================================\n"
			<< "-------------------------  M4 SOLVER RUNNING  --------------------------\n"
			<< "........................................................................\n";
		std::cout
			<< "Current time: " << theWatch().getDayTime()
			<< "    Date: " << theWatch().getDate() << "\n";
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
}