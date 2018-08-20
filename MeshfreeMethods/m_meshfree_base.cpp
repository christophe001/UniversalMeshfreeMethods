/*! \file m_meshfree_base.cpp */

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

#include "m_meshfree_base.h"
#include "u_timer.h"
#include <windows.h>
#include <fstream>
#include <cmath>

constexpr auto theWatch = msl::Timer::instance;

namespace msl {
	MeshfreeBase::MeshfreeBase(double h2dp, bool const_ts, bool enable_rx, double rc) 
		: horizon_ratio_(h2dp), constant_timestep_(const_ts), relaxation_(enable_rx), rc_(rc) {
		class_name_ = "MeshfreeBase";
		sv_step_ = 1;
		sv_count_ = 0;
	}

	void MeshfreeBase::createEnsemble(ObjInfo obj, DomainConfig domain_config) {
		dp_ = obj.dp;
		horizon_ = dp_ * horizon_ratio_;
		density_ = obj.density;
		std::cout
			<< "*************************  CREATING ENSEMBLE  **************************\n"
			<< "creating..\n";
		creator_ = std::make_shared<EnsembleCreator>(obj.shape, obj.offset);
		std::vector<std::string> sas;
		std::vector<std::string> vas;
		std::vector<std::string> tas;
		AttributesGen gen(obj);
		auto infos = gen.getAttrs();
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
		creator_->setDims(obj.dims);
		creator_->setOrientation(obj.orientation, obj.theta);
		creator_->setDpDensity(dp_, density_);
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
		ensemble_->calcBox();
		pos_ = ensemble_->getPos();
		vel_ = ensemble_->getVel();
		acc_ = ensemble_->getAcc();
		np_ = ensemble_->getNp();
		id_ = ensemble_->getId();
		dict_ = ensemble_->getDict();
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
		configSolver(obj.solver);
		setInitVel(obj.velocity);
	}

	void MeshfreeBase::printLabel() const {
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

	void MeshfreeBase::printInfo() const {
		std::cout
			<< "***************************  FILE LOCATION  ****************************\n"
			<< "Target folder: " << folder_
			<< "\t File name: " << filename_ << std::endl
			<< "\n";
		std::cout
			<< "************************  DOMAIN CONFIGURATION  ************************\n"
			<< domain_config_ << std::endl << "\n"
			<< "*************************  ENSEMBLE INFOMATION  ************************\n"
			<< "Total particles: " << np_ << std::endl
			<< "Inter-particle distance: " << dp_ << std::endl
			<< "Pos min: " << (ensemble_->getPosMin()).format(PureFmt) << std::endl
			<< "Pos max: " << (ensemble_->getPosMax()).format(PureFmt) << std::endl
			<< "\n";
		std::cout
			<< "************************  COMPUTATION DETAILS  *************************\n";
		if (l_compute_ != nullptr) {
			std::cout << "Lagrangian Solver: " << l_compute_->getClassName() << std::endl;
			std::cout << "Lagrangian compute params: " << l_compute_->format();
		}
		if (e_compute_ != nullptr) {
			std::cout << "Eulerian solver: " << e_compute_->getClassName() << std::endl;
			std::cout << "Lagrangian compute params: " << e_compute_->format();
		}
		std::cout << std::endl;
	}

	void MeshfreeBase::adaptiveRelaxation() {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
			for (int i = 0; i < np_; i++)
				acc_[i] -= rc_ * vel_[i];
	}

	void MeshfreeBase::setInitVel(Vec3d vel) {
		init_vel_ = vel;
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++) {
			vel_[i] = vel;
		}
	}

	void MeshfreeBase::addPosForce(std::shared_ptr<Shape> shape, Vec3d force) {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++) {
			if (shape->isWithin(pos_[i])) {
				acc_[i] += force;
			}
		}
	}

	void MeshfreeBase::configPosConstraint(std::shared_ptr<Shape> shape) {
		confined_regions_.push_back(shape);
	}

	void MeshfreeBase::configPosConstraintComplement(std::shared_ptr<Shape> shape) {
		confined_regions_complement_.push_back(shape);
	}

	void MeshfreeBase::addInitialPosForce(std::shared_ptr<Shape> shape, Vec3d force) {
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

	void MeshfreeBase::configExternalPosForce(std::shared_ptr<Shape> shape, Vec3d force) {
		force_regions_.push_back(shape);
		forces_.push_back(force);
	}

	void MeshfreeBase::configExternalInitPosForce(std::shared_ptr<Shape> shape, Vec3d force) {
		force_init_regions_.push_back(shape);
		forces_init_.push_back(force);
	}

	void MeshfreeBase::initForces() {
		for (int i = 0; i < force_regions_.size(); i++) {
			addPosForce(force_regions_[i], forces_[i]);
		}
		for (int i = 0; i < force_init_regions_.size(); i++) {
			addInitialPosForce(force_init_regions_[i], forces_init_[i]);
		}
	}

	void MeshfreeBase::setDt(double dt) {
		if (!constant_timestep_) {
			dt_ = min(dt, dt_max_);
		}
	}

	void MeshfreeBase::printMemoryInfo() const {
		double nb_mem = (nbh_ != nullptr) ? nbh_->getMemoryInMB() : 0.0;
		double pb_mem = (pbh_ != nullptr) ? pbh_->getMemoryInMB() : 0.0;
		std::cout
			<< "***************************  MEMORY SUMMARY  ***************************\n"
			<< "Allocated memory for particle: " << ensemble_->getMemoryInMB() << " MB" << std::endl
			<< "Allocated memory for sorting:  " << sorted_->getMemoryInMB() << " MB" << std::endl
			<< "Allocated memory for bonds: " << nb_mem + pb_mem << " MB" << std::endl
			<< "Total allocated memory: " << ensemble_->getMemoryInMB() + sorted_->getMemoryInMB() + nb_mem + pb_mem 
			<< " MB" << std::endl
			<< "\n";
	}

	void MeshfreeBase::configSolver(std::string lcompute, std::string ecompute, 
		std::vector<std::string> attrs) {
		if (!lcompute.empty()) {
			std::shared_ptr<ComputeNeighbor> cpn = std::make_shared<ComputeNeighbor>(sorted_);
			//std::cout << "dp: " << dp_ << " horizon: " << horizon_ << std::endl;
			cpn->initialize(dp_, horizon_, attrs);
			cpn->compute();
			nbh_ = cpn->getNeighborhoodData();
			pbh_ = cpn->getPeriNeighborData();
		}
		if (lcompute == "kl")
			configLagrangianKL();
		if (lcompute == "rigid")
			configRigidBody();
		std::cout << "solver configuration succeed!\n\n";
	}

	void MeshfreeBase::configLagrangianKL() {
		l_compute_ = std::make_shared<StateKLModel>(sorted_, nbh_, pbh_);
	}

	void MeshfreeBase::configRigidBody() {
		l_compute_ = std::make_shared<RigidBody>(sorted_);
	}

	void MeshfreeBase::verletUpdate() {
		if (!constant_timestep_)
			configDt();
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
		for (int i = 0; i < np_; i++)
			pos_[i] += vel_[i] * dt_ + 0.5 * dt_ * dt_ * stack_[i];
		if (!confined_regions_.empty()) {
			Vec3d * ipos = ensemble_->getVectorAttrPtr("initial_position")->getAttr();
			for (auto shape : confined_regions_) {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
				for (int i = 0; i < np_; i++) {
					if (shape->isWithin(ipos[i])) {
						pos_[dict_[i]] = ipos[i];
					}
				}
			}
		}
		if (!confined_regions_complement_.empty()) {
			Vec3d * ipos = ensemble_->getVectorAttrPtr("initial_position")->getAttr();
			for (auto shape : confined_regions_complement_) {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
				for (int i = 0; i < np_; i++) {
					if (!shape->isWithin(ipos[i])) {
						pos_[dict_[i]] = ipos[i];
					}
				}
			}
		}
		initForces();
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

#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
		for (int i = 0; i < np_; i++) {
			stack_[i] = acc_[i];
			acc_[i] = Vec3d::Zero();
		}
		//ensemble_->calcBox();
	}

	void MeshfreeBase::saveVtk() {
		std::string inner_folder = folder_ + "\\" + filename_;
		std::string vtkname = filename_ + "_" + std::to_string(sv_count_) + ".vtk";
		//double curr_time = timestep_ * dt_;

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

		file << "VECTORS cell_id double\n";
		for (int i = 0; i < np_; i++) 
			file << domain_config_.getCellNum(pos_[i]).format(PureFmt) << "\n";
			
		file << "VECTORS Vel double\n";
		for (int i = 0; i < np_; i++) {
			file << vel_[i].format(PureFmt) << "\n";
		}

		file << "VECTORS Acc double\n";
		for (int i = 0; i < np_; i++) {
			file << stack_[i].format(PureFmt) << "\n";
		}

		for (auto&& s : vector_saves_) {
			auto va = ensemble_->getVectorAttrPtr(s);
			file << "VECTORS " << va->getInfo() << " double\n";
			for (int i = 0; i < np_; i++)
				file << va->getAttr()[i].format(PureFmt) << "\n";
		}
		sv_count_++;
		file.close();
	}

	void MeshfreeBase::configIO(std::string folder, std::string filename) {
		folder_ = folder;
		filename_ = filename;
		std::string inner_folder = folder_ + "\\" + filename_;
		CreateDirectory(std::wstring(folder_.begin(), folder_.end()).c_str(), NULL);
		CreateDirectory(std::wstring(inner_folder.begin(), inner_folder.end()).c_str(), NULL);
	}

	void MeshfreeBase::configRun(double dt, double total_time, int sv_step) {
		sv_step_ = sv_step;
		dt_max_ = dt_ = dt;
		T_ = total_time;
		t_intv_ = sv_step_ * dt_;
		t_cur_ = 0.0;
		timestep_ = 0;
		sim_time_ = constant_timestep_ ? (int(T_ / dt_) + 1) * dt_ : (int(T_ / t_intv_) + 1) * t_intv_;
	}

	void MeshfreeBase::run() {
		if (l_compute_ == nullptr && e_compute_ == nullptr) return;
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
			for (; timestep_ <= int(T_/dt_); timestep_++) {
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
			for (; timestep_ <= int(T_ / dt_); timestep_++) {
				saveVtk();
				std::cout << " SaveVtk: " << sv_count_ << "\tTimestep: " << timestep_
					<< "\t" << std::endl;
				while(t_cur_ < t_intv_) {
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

	std::vector<AttributeInfo> AttributesGen::getAttrs() {
		if (obj_.solver == "rigid")
			return getRigid();
		if (obj_.solver == "kl")
			return getStateKL();

	}

	std::vector<AttributeInfo> AttributesGen::getRigid() {
		return std::vector<AttributeInfo>{AttributeInfo("vector", "initial_position", true, "pos")};
	}

	std::vector<AttributeInfo> AttributesGen::getStateKL(){
		auto attrs = std::vector<AttributeInfo>();
		typedef AttributeInfo MyAttr;
		attrs.push_back(MyAttr("tensor", "deformation", false, "identity"));
		attrs.push_back(MyAttr("tensor", "deformation_last", false, "identity"));
		attrs.push_back(MyAttr("tensor", "shape_tensor", false, "zero"));
		attrs.push_back(MyAttr("tensor", "deformation_dot", false, "zero"));
		attrs.push_back(MyAttr("tensor", "hencky", false, "zero"));
		//attrs.push_back(MyAttr("tensor", "rotation", false, "zero"));
		//attrs.push_back(MyAttr("tensor", "left_stretch", false, "zero"));
		//attrs.push_back(MyAttr("tensor", "left_stretch", false, "zero"));
		//attrs.push_back(MyAttr("tensor", "left_stretch_dot", false, "zero"));
		//attrs.push_back(MyAttr("tensor", "d", false, "zero"));
		attrs.push_back(MyAttr("tensor", "tau", false, "zero"));
		attrs.push_back(MyAttr("vector", "initial_position", true, "pos"));
		attrs.push_back(MyAttr("scalar", "density", true, "constant", obj_.density));
		attrs.push_back(MyAttr("scalar", "pressure", true, "zero"));
		attrs.push_back(MyAttr("scalar", "damage", true, "zero"));
		attrs.push_back(MyAttr("scalar", "delta_p", false, "zero"));
		attrs.push_back(MyAttr("scalar", "theta_e", true, "zero"));
		attrs.push_back(MyAttr("scalar", "theta_p", true, "zero"));
		return attrs;
	}
}