/*! \file m_state_peridynamics.cpp */

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

#include "m_state_peridynamics.h"
#include <unsupported/Eigen/MatrixFunctions>
#ifdef _DEBUG_
#include <iostream>
#include <random>
#include "u_timer.h"
#endif // _DEBUG_
namespace msl {
	StateBasedPD::StateBasedPD(std::shared_ptr<SortEnsemble> sorted,
		std::shared_ptr<NeighborhoodData> nbh, std::shared_ptr<PeriNeighborData> pbh, 
		double h2dp, double dt, bool damage) 
		: LagrangianCompute::LagrangianCompute(sorted) {
		//LagrangianCompute::computeBond();
		nbh_ = nbh;
		pnbh_ = pbh;
		class_name_ = "StateBasedPD";
		np_ = ensemble_ptr_->getNp();
		damage_enabled_ = damage;
		force_const_ = 1.5 * pow(10, 9);
		shape_tensor_		= ensemble_ptr_->getTensorAttrPtr("shape_tensor")->getAttr();
		deformation_		= ensemble_ptr_->getTensorAttrPtr("deformation")->getAttr();
		tau_				= ensemble_ptr_->getTensorAttrPtr("tau")->getAttr();
		if (ensemble_ptr_->hasTensorAttribute("deformation_dot"))
			deformation_dot_ = ensemble_ptr_->getTensorAttrPtr("deformation_dot")->getAttr();
		if (ensemble_ptr_->hasTensorAttribute("deformation_last"))
			deformation_last_ = ensemble_ptr_->getTensorAttrPtr("deformation_last")->getAttr();
		if (ensemble_ptr_->hasTensorAttribute("rotation"))
			rotation_ = ensemble_ptr_->getTensorAttrPtr("rotation")->getAttr();
		if (ensemble_ptr_->hasTensorAttribute("left_stretch"))
			left_stretch_ = ensemble_ptr_->getTensorAttrPtr("left_stretch")->getAttr();
		if (ensemble_ptr_->hasTensorAttribute("left_stretch_dot"))
			left_stretch_dot_ = ensemble_ptr_->getTensorAttrPtr("left_stretch_dot")->getAttr();
		if (ensemble_ptr_->hasTensorAttribute("d"))
			d_ = ensemble_ptr_->getTensorAttrPtr("d")->getAttr();
		init_pos_	= ensemble_ptr_->getVectorAttrPtr("initial_position")->getAttr();
		density_	= ensemble_ptr_->getScalarAttrPtr("density")->getAttr();
		if (ensemble_ptr_->hasScalarAttribute("lame_lambda"))
			lambda_ = ensemble_ptr_->getScalarAttrPtr("lame_lambda")->getAttr();
		if (ensemble_ptr_->hasScalarAttribute("lame_mu"))
			mu_ = ensemble_ptr_->getScalarAttrPtr("lame_mu")->getAttr();
		if (damage_enabled_) {
			if (ensemble_ptr_->hasScalarAttribute("damage1"))
				damage_ = ensemble_ptr_->getScalarAttrPtr("damage1")->getAttr();
			else throwException("Constructor", "No damage1 variable found!");
			if (ensemble_ptr_->hasScalarAttribute("damage"))
				damage_all_ = ensemble_ptr_->getScalarAttrPtr("damage")->getAttr();
			else throwException("Constructor", "No damage variable found!");
		}
		else damage_ = nullptr;

		nbl_ = nbh_->getNeighborhoodList();
		
		i_epsilon_ = 0;
		i_p_ = 1;
		shape_calc_ = false;
		dt_ = dt;
		dp_ = sorted->getEnsemble()->getDp();
		rho_ = sorted->getEnsemble()->getDensity();
		horizon_ = dp_ * h2dp;


	}

	StateBasedPD::StateBasedPD(std::shared_ptr<ComputeNeighbor> cpn)
		: LagrangianCompute::LagrangianCompute(cpn) {
		shape_tensor_ = ensemble_ptr_->getTensorAttrPtr("shape_tensor")->getAttr();
		deformation_ = ensemble_ptr_->getTensorAttrPtr("deformation")->getAttr();
		deformation_dot_ = ensemble_ptr_->getTensorAttrPtr("deformation_dot")->getAttr();
		rotation_ = ensemble_ptr_->getTensorAttrPtr("rotation")->getAttr();
		left_stretch_ = ensemble_ptr_->getTensorAttrPtr("left_stretch")->getAttr();
		left_stretch_dot_ = ensemble_ptr_->getTensorAttrPtr("left_stretch_dot")->getAttr();
		d_ = ensemble_ptr_->getTensorAttrPtr("d")->getAttr();
		tau_ = ensemble_ptr_->getTensorAttrPtr("tau")->getAttr();
		init_pos_ = ensemble_ptr_->getVectorAttrPtr("initial_position")->getAttr();
		density_ = ensemble_ptr_->getScalarAttrPtr("density")->getAttr();
		lambda_ = ensemble_ptr_->getScalarAttrPtr("lame_lambda")->getAttr();
		mu_ = ensemble_ptr_->getScalarAttrPtr("lame_mu")->getAttr();
		nbl_ = nbh_->getNeighborhoodList();
		np_ = ensemble_ptr_->getNp();
		i_epsilon_ = 0;
		i_p_ = 2;
		shape_calc_ = false;
		class_name_ = "StateBasedPD";
	}

	void StateBasedPD::configModelParams(const std::vector<std::vector<double>>& params) {
		if (params.size() == 0)
			throwException("configModelParams", "no params found");
		if (params[0].size() < 4)
			throwException("configModelParams", "state based PD not enough params");
		configParams(params[0][0], params[0][1], params[0][2], params[0][3]);
	}

	std::string StateBasedPD::format() const {
		std::string s;
		s += ("  density: " + std::to_string(rho_) + "\n");
		s += ("  horizon: " + std::to_string(horizon_) + "\n");
		s += ("  time_step: " + std::to_string(dt_) + "\n");
		s += ("  dp: " + std::to_string(dp_) + "\n");
		s += ("  influence function parameters: {epsilon: "
			+ std::to_string(i_epsilon_) + " }, { p: " + std::to_string(i_p_) + " }\n");
		return s;
	}

	void StateBasedPD::configParams(double horizon, double dp, double dt, double rho) {
		dp_ = dp;
		epsilon_ = dp_;
		horizon_ = horizon;
		dt_ = dt;
		rho_ = rho;
		dv_ = sorted_ptr_->getDomainConfig()->sim2d() ? dp_ * dp_ : dp_ * dp_ * dp_;
	}
	
	double StateBasedPD::volumeCorrector(double d) {
		return d < horizon_ - dp_ / 2.0 ? 1 : (horizon_ + dp_ / 2.0 - d) / dp_;
	}

	void StateBasedPD::addShapeTensor(int i, long it) {
		if (nbh_->getBondDamage()[it] < 1.0) {
			int j = nbl_[it];
			Vec3d xi = init_pos_[j] - init_pos_[i];
			double d = xi.norm();
			shape_tensor_[i] += influence(d) * correction(i, j) * xi * xi.transpose() * volumeCorrector(d);
		}
	}

	void StateBasedPD::addDeformation(int i, long it) {
		if (nbh_->getBondDamage()[it] < 1.0) {
			int j = nbl_[it];
			Vec3d xi = init_pos_[j] - init_pos_[i];
			Vec3d eta = pos_[dict_[j]] - pos_[dict_[i]];
			double d = xi.norm();
			deformation_[i] += influence(d) * correction(i, j) * eta * xi.transpose() * volumeCorrector(d);
		}
	}

	void StateBasedPD::addDeformationGradient(int i, long it) {
		if (nbh_->getBondDamage()[it] < 1.0) {
			int j = nbl_[it];
			Vec3d xi = init_pos_[j] - init_pos_[i];
			Vec3d eta_dot = vel_[dict_[j]] - vel_[dict_[i]];
			double d = xi.norm();
			deformation_dot_[i] += influence(d) * correction(i, j) * eta_dot * xi.transpose() * volumeCorrector(d);
		}
	}

	void StateBasedPD::addForceState(int i, long it) {
		if (nbh_->getBondDamage()[it] < 1.0) {
			int j = nbl_[it];
			Vec3d xi = init_pos_[j] - init_pos_[i];
			double d = xi.norm();
			Mat3d sigma_i = rotation_[i] * tau_[i] * rotation_[i].transpose();
			Mat3d sigma_j = rotation_[j] * tau_[j] * rotation_[j].transpose();
			Mat3d P_i = deformation_[i].determinant() * sigma_i * deformation_[i].inverse().transpose();
			Mat3d P_j = deformation_[j].determinant() * sigma_j * deformation_[j].inverse().transpose();
			Vec3d force_i = influence(d) * correction(i, j) * (P_i * shape_tensor_[i] + P_j * shape_tensor_[j]) * xi;
			acc_[dict_[i]] += force_i * volumeCorrector(d) * pow(10, 9) / rho_;
		}
	}

	void StateBasedPD::addForceStateNoRot(int i, long it) {
		if (nbh_->getBondDamage()[it] < 1.0) {
			int j = nbl_[it];
			if (damage_[i] < 0.99) {
				Vec3d xi = init_pos_[j] - init_pos_[i];
				double d = xi.norm();
				Mat3d P_i = deformation_[i].determinant() * tau_[i] * deformation_[i].inverse().transpose();
				Mat3d P_j = deformation_[j].determinant() * tau_[j] * deformation_[j].inverse().transpose();
				Vec3d force_i = influence(d) * correction(i, j) * (P_i * shape_tensor_[i] + P_j * shape_tensor_[j]) * xi;
				acc_[dict_[i]] += force_i * volumeCorrector(d) * pow(10, 9) / rho_;
			}
			
			if ((pos_[dict_[j]] - pos_[dict_[i]]).norm() < epsilon_) {
				Vec3d eta = pos_[dict_[j]] - pos_[dict_[i]];
				acc_[dict_[i]] -= force_const_ * dv_ * eta / (eta.norm() * ensemble_ptr_->getMass())
					* log(epsilon_ / eta.norm());
			}
		}
	}

	void StateBasedPD::computeShapeTensor() {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
			for (int i = 0; i < np_; i++)
				shape_tensor_[i] = Mat3d::Zero();
			compute(static_cast<LagrangianCompute::fpn>(&StateBasedPD::addShapeTensor));
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
			for (int i = 0; i < np_; i++) {
				if (shape_tensor_[i].determinant() == 0)
					shape_tensor_[i] = Mat3d::Identity();
				else { 
					Mat3d temp = shape_tensor_[i].inverse(); 
					shape_tensor_[i] = temp;
				}
			}
#ifdef _DEBUG_
			std::uniform_int_distribution<> dis(0, np_ - 1);
			std::default_random_engine re;
			re.seed(std::chrono::system_clock::now().time_since_epoch().count());
			for (int i = 0; i < 5; i++) {
				int id = dis(re);
				std::cout << "id: " << id << std::endl
					<< shape_tensor_[id].format(CleanFmt) << "\n\n";
			}
#endif // !_DEBUG_

	}

	void StateBasedPD::computeDeformation() {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
		for (int i = 0; i < np_; i++)
			deformation_[i] = Mat3d::Zero();

		compute(static_cast<LagrangianCompute::fpn>(&StateBasedPD::addDeformation));
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
		for (int i = 0; i < np_; i++) {
			if (deformation_[i].isZero())
				deformation_[i] = deformation_last_[i];
			else
				deformation_[i] = deformation_[i] * shape_tensor_[i];			
			//
			if (deformation_[i].determinant() == 0.0)
				deformation_[i] += 0.01 * Mat3d::Identity();
			//
		}
#ifdef _DEBUG_
		std::uniform_int_distribution<> dis(0, np_ - 1);
		std::default_random_engine re;
		re.seed(std::chrono::system_clock::now().time_since_epoch().count());
		for (int i = 0; i < 5; i++) {
			int id = dis(re);
			std::cout << "id: " << id << std::endl
				<< deformation_[id].format(CleanFmt) << "\n\n";
		}
#endif // !_DEBUG_
	}

	void StateBasedPD::computeDeformationDot() {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
		for (int i = 0; i < np_; i++)
			deformation_dot_[i] = Mat3d::Zero();
		compute(static_cast<LagrangianCompute::fpn>(&StateBasedPD::addDeformationGradient));
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
		for (int i = 0; i < np_; i++)
			deformation_dot_[i] = deformation_dot_[i] * shape_tensor_[i];
	}


	void StateBasedPD::updateLeftStretch() {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
		for (int i = 0; i < np_; i++)
			left_stretch_[i] += dt_ * left_stretch_dot_[i];
	}

	void StateBasedPD::updateDensity() {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
		for (int i = 0; i < np_; i++)
			density_[i] = rho_ / deformation_[i].determinant();
	}

	void StateBasedPD::computeForces() {
		computeShapeTensor();
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++)
			deformation_last_[i] = deformation_[i];
		computeDeformation();
		computeDeformationDot();
		updateDensity();
		computeRotation();
		computeTau();
		compute(static_cast<LagrangianCompute::fpn>(&StateBasedPD::addForceState));
	}


	void StateBasedPD::computeRotation() {
		updateLeftStretch();
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
		for (int i = 0; i < np_; i++) {
			Mat3d L = deformation_dot_[i] * deformation_[i].inverse();
			Mat3d V = left_stretch_[i];
			Mat3d D = 0.5 * (L + L.transpose());
			Mat3d W = 0.5 * (L - L.transpose());
			Mat3d T = D * W;
			Vec3d w = -0.5 * Vec3d(W(1, 2) - W(2, 1), W(2, 0) - W(0, 2), W(0, 1) - W(1, 0));
			Vec3d z(T(1, 2) - T(2, 1), T(2, 0) - T(0, 2), T(0, 1) - T(1, 0));
			Mat3d tv = Mat3d::Identity() * V.trace() - V;
			Vec3d omega = w;
			if (tv.determinant() != 0)
				omega += tv.inverse() * z;
			Mat3d Omega;
			Omega << 0.0, -omega(2), omega(1),
				omega(2), 0.0, -omega(0),
				-omega(1), omega(0), 0.0;
			left_stretch_dot_[i] = L * V - V * Omega;
			double Omega_square = omega(0)*omega(0) + omega(1)*omega(1) + omega(2)*omega(2);
			double Omega_square_root = sqrt(Omega_square);
			Mat3d d_rotation = Mat3d::Identity();
			if (Omega_square != 0) {
				d_rotation += sin(dt_*Omega_square_root) / Omega_square_root*Omega
					- (1 - cos(dt_*Omega_square_root)) / Omega_square * Omega * Omega;
			}
			rotation_[i] = d_rotation * rotation_[i];
			d_[i] = rotation_[i].transpose() * D * rotation_[i];
		}
	}

}