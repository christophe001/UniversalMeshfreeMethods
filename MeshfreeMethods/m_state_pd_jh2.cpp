/*! \file m_state_pd_jh2.cpp */

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

#include "m_state_pd_jh2.h"
#include <cmath>
#include <algorithm>

namespace msl {
	StateBasedPDJH2::StateBasedPDJH2(std::shared_ptr<SortEnsemble> sorted, 
		std::shared_ptr<NeighborhoodData> nbh, std::shared_ptr<PeriNeighborData> pbh)
		: StateBasedPD(sorted, nbh, pbh) {
		damage_ = ensemble_ptr_->getScalarAttrPtr("damage")->getAttr();
		pressure_ = ensemble_ptr_->getScalarAttrPtr("pressure")->getAttr();
		delta_p_ = ensemble_ptr_->getScalarAttrPtr("delta_p")->getAttr();
		class_name_ = "StateBasedPDJH2";
	}

	StateBasedPDJH2::StateBasedPDJH2(std::shared_ptr<ComputeNeighbor> cpn)
		: StateBasedPD(cpn) {
		damage_ = ensemble_ptr_->getScalarAttrPtr("damage")->getAttr();
		pressure_ = ensemble_ptr_->getScalarAttrPtr("pressure")->getAttr();
		delta_p_ = ensemble_ptr_->getScalarAttrPtr("delta_p")->getAttr();
		class_name_ = "StateBasedPDJH2";
	}


	void StateBasedPDJH2::computeTau() {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP
		for (int i = 0; i < np_; i++) {
			// Compute pressure
			double mu = density_[i] / rho_ - 1;
			if (mu >= 0) {
				pressure_[i] = delta_p_[i] + K1_ * mu + K2_ * mu * mu + K3_ * pow(mu, 3);
			}
			else {
				pressure_[i] = delta_p_[i] + K1_ * mu;
			}
			// Strain update
			Mat3d delta_e = d_[i] * dt_;
			// Compute trial stress
			Mat3d tau_tr = tau_[i] + lambda_[i] * delta_e.trace() * Mat3d::Identity() + 2 * mu_[i] * delta_e;
			// Deviatoric trial stress
			Mat3d S_tr = tau_tr - 1.0 / 3 * tau_tr.trace()*Mat3d::Identity();
			// Equivalent stress
			double S_vm = sqrt(1.5 * (S_tr.transpose() * S_tr).trace());
			// Effective strain rate
			double d_vm = sqrt(1.5 * (d_[i].transpose() * d_[i]).trace());
			// Compute normalized stress
			double p_s = pressure_[i] / p_hel_;
			double T_s = T_ / p_hel_;
			double sigma_i = A_ * pow(p_s + T_s, N_) * (1 + C_ * log(d_vm));
			double sigma_f = B_ * pow(p_s, M_) * (1 + C_ * log(d_vm));
			double sigma_s = sigma_i - damage_[i] * (sigma_i - sigma_f);
			double epsilon_f = D1_ * pow(p_s + T_s, D2_);
			// Check for yielding
			if (sigma_s * sigma_hel_ < S_vm) {
				// Compute Lagrangian multiplier
				Mat3d ce = lambda_[i] * d_[i].trace() * Mat3d::Identity() + 2 * mu_[i] * d_[i];
				Mat3d ca = 1.5 / S_vm * (lambda_[i] * S_tr.trace() * Mat3d::Identity() + 2 * mu_[i] * S_tr);
				Mat3d a = 1.5 / S_vm * S_tr;
				double lagmul = (a.cwiseProduct(ce).sum()) / (a.cwiseProduct(ca).sum());
				double sc = pow(S_tr(0, 0) - S_tr(1, 1), 2) + pow(S_tr(1, 1) - S_tr(2, 2), 2) + pow(S_tr(2, 2) - S_tr(0, 0), 2)
					+ 6 * (pow(S_tr(0, 1), 2) + pow(S_tr(0, 2), 2) + pow(S_tr(1, 2), 2));
				// Compute epsilon_p, update damage
				double epsilon_p = lagmul * dt_ * sqrt(0.5 * sc);
				damage_[i] += epsilon_p / epsilon_f;
				damage_[i] = std::min(damage_[i], 1.0);
				double u_1 = pow(sigma_s * sigma_hel_, 2) / (6 * mu_[i]);
				double u_2 = pow(sigma_hel_ * (sigma_i - damage_[i] * (sigma_i - sigma_f)), 2) / (6 * mu_[i]);
				// Update delta_p
				double delta_u = u_1 - u_2;
				delta_p_[i] = -K1_ * mu + sqrt(-K1_ * mu + delta_p_[i] + 2 * K1_ * delta_u);
				double p_t = 1.0 / 3.0 * tau_tr.trace() + delta_p_[i];
				tau_tr = S_tr + p_t * Mat3d::Identity();
			}
			// Update tau
			tau_[i] = tau_tr;
		}
	}

	void StateBasedPDJH2::configModelParams(const std::vector<std::vector<double>>& params) {
		if (params.size() < 5)
			throwException("configModelParams", "not enough params");
		StateBasedPD::configModelParams(params);
		if (params[1].size() < 5)
			throwException("configParamsA", "not enough params");
		if (params[2].size() < 2)
			throwException("configParamsD", "not enough params");
		if (params[3].size() < 2)
			throwException("configParamsHel", "not enough params");
		if (params[4].size() < 3)
			throwException("configParamsK", "not enough params");
		configParamsA(params[1][0], params[1][1], params[1][2], params[1][3], params[1][4]);
		configParamsD(params[2][0], params[2][1]);
		configParamsHel(params[3][0], params[3][1]);
		configParamsK(params[4][0], params[4][1], params[4][2]);
	}

	std::string StateBasedPDJH2::format() const {
		std::string s = StateBasedPD::format();
		s += "JH2 model parameters:\n";
		s += "  Constants:\n";
		s += "  A: " + std::to_string(A_) + "  B: " + std::to_string(B_) + "  C: " + std::to_string(C_)
			+ "  M: " + std::to_string(M_) + "  N: " + std::to_string(N_) + "\n";
		s += "  HEL params:\n";
		s += "  sigma_HEL: " + std::to_string(sigma_hel_) + "  p_HEL: " + std::to_string(p_hel_) + "\n";
		s += "  Damage constants:\n";
		s += "  D_1: " + std::to_string(D1_) + "  D_2: " + std::to_string(D2_) + "\n";
		s += "  Equation of state:\n";
		s += "  K_1: " + std::to_string(K1_) + "  K_2: " + std::to_string(K2_) + "  K_3: " + std::to_string(K3_) + "\n";
		return s;
	}

}