/*! \file m_state_kl_model.cpp */

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

#include "m_state_kl_model.h"
#include <cmath>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
namespace msl {
	//! core
	void StateKLModel::computeTau() {
		/*
		std::cout << "test\n" << " damage: " << v.damage << std::endl;
		std::cout << "theta e: " << v.theta_e << std::endl;
		std::cout << "theta p: " << v.theta_p << std::endl;
		std::cout << "hencky: \n" << v.hencky << std::endl;
		std::cout << "tau: \n" << v.sigma << std::endl;
		*/
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_
		for (int i = 0; i < np_; i++) {
			auto var = dgamma(tau_[i], hencky_[i], deformation_last_[i], 
				deformation_[i], theta_e_[i], theta_p_[i], damage_[i]);
			damage_[i] = var.damage;
			theta_e_[i] = var.theta_e;
			theta_p_[i] = var.theta_p;
			hencky_[i] = var.hencky;
			tau_[i] = var.sigma;
		}
	}

	void StateKLModel::computeForces() {
#ifdef _WITH_OMP_
#pragma omp parallel for schedule(static)
#endif // _WITH_OMP_

		for (int i = 0; i < np_; i++)
			deformation_last_[i] = deformation_[i];
		StateBasedPD::computeShapeTensor();
		StateBasedPD::computeDeformation();
		//StateBasedPD::computeDeformationDot();
		//StateBasedPD::updateDensity();
		StateKLModel::computeTau();
		compute(static_cast<LagrangianCompute::fpn>(&StateBasedPD::addForceStateNoRot));
	}

	void StateKLModel::setDamageParams(double c1, double c2, double d1, double d2) {
		C1_ = c1;
		C2_ = c2;
		D1_ = d1;
		D2_ = d2;
	}

	StateKLModel::StateKLModel(std::shared_ptr<SortEnsemble> sorted,
		std::shared_ptr<NeighborhoodData> nbh, std::shared_ptr<PeriNeighborData> pbh, double dt)
		: StateBasedPD(sorted, nbh, pbh, dt),
		a_{ { 28.89, -86.62, -467.50, -710.22 },
			{ 5.51, 67.47, -59.49, -297.95 },
			{ 116.22, 3723.38, 18368.38, 21463.56 }}, 
		b_ {{-29.26, 359.97, 1592.37, 2351.28},
			{-104.25, -2558.43, -10752.23, -9194.98}} {
		damage_ = ensemble_ptr_->getScalarAttrPtr("damage")->getAttr();
		pressure_ = ensemble_ptr_->getScalarAttrPtr("pressure")->getAttr();
		theta_e_ = ensemble_ptr_->getScalarAttrPtr("theta_e")->getAttr();
		theta_p_ = ensemble_ptr_->getScalarAttrPtr("theta_p")->getAttr();
		hencky_ = ensemble_ptr_->getTensorAttrPtr("hencky")->getAttr();
		klog_ = false;
		//! material parameters for consolidation
		p_0_ = 16.937;
		Lambda_ = -5.226;
		alpha_ = 9.203;

		//! material parameters for critical state
		p_1_ = 12.337;
		p_n_ = -8.5;
		q_n_ = 7.402;
		B_ = 1.168;
		beta_ = 0.5;
		p_i_ = 4.99;

		//! tensile failure pressure
		p_t_ = -10;

		//! parameters for damage
		C1_ = 0.1;
		C2_ = 0.0;
		D1_ = 0.005;
		D2_ = 0.0;
		class_name_ = "StateKLModel";
	}


	//! Yield condition
	double StateKLModel::yield(const Mat3d& sigma, const double& theta_p, const double& damage) const {
		double p = -sigma.trace() / 3.0;
		double q = sqrt(0.5) * dev(sigma).norm();
		double p_c = pc(theta_p);
		double q_f = qf(p_c);
		double q_c = qc(0.5 * p_c + 0.5 * p_t_);
		return q * q + 4.0 * pow((1 - damage) * q_c + damage * q_f, 2.0) * 
			(p - (1.0 - damage) * p_t_) * (p - p_c) / pow(p_c - (1 - damage) * p_t_, 2.0);
	}

	double StateKLModel::checkYield(const Mat3d& h_tr, const Mat3d& sigma, const double& J, const double& theta_e,
		const double& theta_p, const double& p, const double& damage, const double& dgamma) const {
		auto theta = actualTheta(dgamma, p, theta_e, theta_p, damage);
		double theta_ne = theta.first, theta_np = theta.second;
		double damage_n = actualDamage(sigma, dgamma, damage);
		auto sig_h = actualSigmaHencky(h_tr, J, dgamma, theta_ne, theta_np, p, damage_n);
		return yield(sig_h.first, theta_np, damage_n);
	}


	double StateKLModel::pc(const double& theta_p) const {
		return p_0_ + Lambda_ / alpha_ * (1 - exp(-alpha_ * theta_p));
	}

	double StateKLModel::qc(const double& pm) const {
		if (pm <= p_i_)
			return (p_1_ - pm) / (p_1_ - p_n_) * q_n_;
		else
			return B_ * pow(pm, beta_);
	}

	double StateKLModel::qf(const double& pc) const{
		return C1_ * pow(0.5 * pc, C2_);
	}

	//! Elastic response
	double StateKLModel::shear(const double& theta_e, const double& theta_p) const {
		double a0 = a_[0][0] + a_[0][1] * theta_p + a_[0][2] * pow(theta_p, 2.0) + a_[0][3] * pow(theta_p, 3.0);
		double a1 = a_[1][0] + a_[1][1] * theta_p + a_[1][2] * pow(theta_p, 2.0) + a_[1][3] * pow(theta_p, 3.0);
		double a2 = a_[2][0] + a_[2][1] * theta_p + a_[2][2] * pow(theta_p, 2.0) + a_[2][3] * pow(theta_p, 3.0);
		return a0 + a1 * theta_e + a2 * theta_e * theta_e;
	}

	double StateKLModel::dshear(const double& theta_e, const double& theta_p) const {
		double a1 = a_[1][0] + a_[1][1] * theta_p + a_[1][2] * pow(theta_p, 2.0) + a_[1][3] * pow(theta_p, 3.0);
		double a2 = a_[2][0] + a_[2][1] * theta_p + a_[2][2] * pow(theta_p, 2.0) + a_[2][3] * pow(theta_p, 3.0);
		return a1 + 2.0 * a2 * theta_e;
	}

	double StateKLModel::Jp(const double& theta_e, const double& theta_p) const {
		double b0 = b_[0][0] + b_[0][1] * theta_p + b_[0][2] * pow(theta_p, 2.0) + b_[0][3] * pow(theta_p, 3.0);
		double b1 = b_[1][0] + b_[1][1] * theta_p + b_[1][2] * pow(theta_p, 2.0) + b_[1][3] * pow(theta_p, 3.0);
		return b0 * theta_e + b1 * pow(theta_e, 3.0);
	}

	double StateKLModel::ddamage(const double& damage, const Mat3d& sigma) const {
		if (damage == 1.0)
			return 0.0;
		double p = -sigma.trace() / 3.0;
		double epf = D1_ * pow(p - p_t_, D2_);
		return sqrt(2.0/3.0) * dev(sigma).norm() / epf;
	}

	//! Function
	double StateKLModel::zeta(const double& p, const double& theta_p, const double& damage) const {
		double p_c = pc(theta_p);
		double q_f = qf(p_c);
		double q_c = qc(0.5 * p_c + 0.5 * p_t_);
		return -4.0 * pow((1 - damage) * q_c + damage * q_f, 2.0) *
			(2.0 * p - (p_c + (1 - damage) * p_t_)) / pow(p_c - (1 - damage) * p_t_, 2.0);
	}

	Mat3d StateKLModel::dev(const Mat3d& m) const {
		return m - 1.0 / 3.0 * m.trace() * Mat3d::Identity();
	}

	double StateKLModel::cbar(const double& b_s, const double& b_t) const {
		return (b_t + b_s) / (b_t - b_s) + 2.0 / (log(b_s) - log(b_t));
	}

	double StateKLModel::ctilde(const double& b_s, const double& b_t) const {
		return 2.0 * sqrt(b_s * b_t) / (b_t - b_s) + 2.0 / (log(b_s) - log(b_t));
	}

	//! Return actual updates
	std::pair<double, double> StateKLModel::actualTheta(const double& dgamma,
		const double& p, const double&theta_e_tr, const double& theta_p, const double& damage) const {
		double c = dgamma * zeta(p, theta_p, damage);
		return std::make_pair(theta_e_tr - c, theta_p + c);
	}

	std::pair<Mat3d, Mat3d> StateKLModel::actualSigmaHencky(const Mat3d& h_tr, const double& J,
		const double& dgamma, const double&theta_e, const double& theta_p, const double& p, const double& damage) const {
		double d = 1.0 - heaviside(-p) * damage;
		double c1 = 2.0 * d * shear(theta_e, theta_p);
		double div = J + dgamma * c1;
		Mat3d sigma = c1 / div * dev(h_tr);
		double pn = -d / J * (dshear(theta_e, theta_p) * pow(J / div, 2.0) * dev(h_tr).squaredNorm() - Jp(theta_e, theta_p));
		Mat3d hencky = h_tr - dgamma * sigma - 1.0 / 3.0 * dgamma * zeta(p, theta_p, damage) * Mat3d::Identity();
		sigma -= pn * Mat3d::Identity();
		return std::make_pair(sigma, hencky);
	}

	double StateKLModel::actualDamage(const Mat3d& sigma, const double& dgamma, const double& damage) const {
		if (damage == 1.0) return 1.0;
		else {
			return std::min(1.0, damage + dgamma * ddamage(damage, sigma));
		}
	}

	//! Trial functions
	Mat3d StateKLModel::trialHencky(const Mat3d& f_last, const Mat3d& f, const Mat3d& b_e) const {
		Mat3d fb = f * f_last.inverse();
		return 0.5 * (fb * b_e * fb.transpose()).log();
	}

	Mat3d StateKLModel::trialSigma(const Mat3d& hencky, const double& J, const double& theta_e,
		const double& theta_p, const double& p, const double& damage) const {
		double a = (1 - heaviside(-p) * damage) / J;
		Mat3d dev_h = dev(hencky);
		Mat3d dev = a * (2.0 * shear(theta_e, theta_p) * dev_h);
		double p_trial = -a * (dshear(theta_e, theta_p) * dev_h.squaredNorm() - Jp(theta_e, theta_p));
		return -p_trial * Mat3d::Identity() + dev;
	}

	//! plastic multiplier delta gamma calculation
	StateKLModel::Vars StateKLModel::dgamma(const Mat3d& sigma, const Mat3d& hencky, const Mat3d& f_last, 
		const Mat3d& f, const double& theta_e, const double& theta_p, const double& damage) const {
		//! trial state
		double J = f.determinant();
		Mat3d b = (2.0 * hencky).exp();
		double p = -sigma.trace() / 3.0;
		Mat3d h_tr = trialHencky(f_last, f, b);
		double te_tr = h_tr.trace();
		Mat3d sig_tr = trialSigma(h_tr, J, te_tr, theta_p, p, damage);
		//! check for yield
		if (yield(sig_tr, theta_p, damage) <= epsilon) {
			return Vars(0.0, te_tr, theta_p, sig_tr, h_tr);
		}
		else {
			double dg = 0.0;
			double der = 0.0;
			while (true) {
				int cnt = 0;
				double c1 = checkYield(h_tr, sigma, J, theta_e, theta_p, p, damage, dg);
				double c2 = checkYield(h_tr, sigma, J, theta_e, theta_p, p, damage, dg + 0.0001);
				if (c1 != c2)
					der = (c2 - c1) / 0.0001;
				else {
					c2 = checkYield(h_tr, sigma, J, theta_e, theta_p, p, damage, dg + 0.00005);
					der = (c2 - c1) / 0.00005;
				}
				if (abs(c1 / der) < pow(10, -6)) {
					break;
				}
				else {
					dg -= c1 / der;
				}
				if (++cnt == 100) {
					std::cout << "dead loops?\n";
					break;
				}
			}
			double d = actualDamage(sigma, dg, damage);
			auto theta_a = actualTheta(dg, p, theta_e, theta_p, damage);
			auto sig_h = actualSigmaHencky(h_tr, J, dg, theta_a.first, theta_a.second, p, d);
			return Vars(d, theta_a.first, theta_a.second, sig_h.first, sig_h.second);
		}
	}

}