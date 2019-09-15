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
#include <Eigen/Eigenvalues> 
#include <random>
#pragma optimize("", off)
using std::max;
using std::min;
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
			try {
				if (damage_[i] < 1.0) {
					StateKLModel::Vars var = dgamma(tau_[i], hencky_[i], deformation_last_[i],
						deformation_[i], theta_e_[i], theta_p_[i], damage1_[i]);
					damage1_[i] = var.damage;
					damage2_[i] = min(1.0, max(damage2_[i],  tensileDamage(i)));
					damage_[i] = 1.0 - (1.0 - damage1_[i]) * (1.0 - damage2_[i]);
					theta_e_[i] = var.theta_e;
					theta_p_[i] = var.theta_p;
					hencky_[i] = var.hencky;
					tau_[i] = var.sigma;
					pressure_[i] = -tau_[i].trace() / 3.0;
				}
			}
			catch (std::exception& e) {
				std::cerr << e.what() << std::endl;
			}
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
		StateBasedPD::updateDensity();
		StateKLModel::computeTau();
		compute(static_cast<LagrangianCompute::fpn>(&StateBasedPD::addForceStateNoRot));

	}

	void StateKLModel::setDamageParams(double c1, double c2, double d1, double d2) {
		C1_ = c1;
		C2_ = c2;
		D1_ = d1;
		D2_ = d2;
	}

	double StateKLModel::tensileDamage(int i) const {
		Mat3d sigma = tau_[i];
		Vec3d ev = sigma.eigenvalues().real();
		double t = ev.maxCoeff();
		return (t - 0.5 * tensile_limit_[i]) / (0.5 * tensile_limit_[i]);
	}

	StateKLModel::StateKLModel(std::shared_ptr<SortEnsemble> sorted, std::shared_ptr<NeighborhoodData> nbh, 
		std::shared_ptr<PeriNeighborData> pbh, double h2dp, double dt)
		: StateBasedPD(sorted, nbh, pbh, h2dp, dt),
		/*
		a_{ { 28.89, -86.62, -467.50, -710.22 },
			{ 5.51, 67.47, -59.49, -297.95 },
			{ 116.22, 3723.38, 18368.38, 21463.56 }},*/
		a_{ {29.85, -45.22861593, 397.29726482, -529.04822999},
			{-70.80208598, 482.7133639, -92.22520802, -2421.98089823},
			{120.13192662, 67.04966709, -5981.958835, 17866.5140235 } },
		/*
		b_ {{-29.26, 359.97, 1592.37, 2351.28},
			{-104.25, -2558.43, -10752.23, -9194.98}},*/
		b_{ {56.60444606, -364.27307119, 1720.75846719, -1817.70208959},
			{-198.05259056, 2786.73610983, -12721.4365674, 14065.20743152},
			{210.04584379, -2935.9322371, 11923.63577417, -2318.61448899} },
		c_ { 3.814, -0.382, 3.16*pow(10.0,-2), -8.03*pow(10,-4), 9.35*pow(10,-6), -4.02*pow(10, -8) }
	{
		damage1_ = ensemble_ptr_->getScalarAttrPtr("damage1")->getAttr();
		damage2_ = ensemble_ptr_->getScalarAttrPtr("damage2")->getAttr();
		damage_ = ensemble_ptr_->getScalarAttrPtr("damage")->getAttr();
		pressure_ = ensemble_ptr_->getScalarAttrPtr("pressure")->getAttr();
		theta_e_ = ensemble_ptr_->getScalarAttrPtr("theta_e")->getAttr();
		theta_p_ = ensemble_ptr_->getScalarAttrPtr("theta_p")->getAttr();
		hencky_ = ensemble_ptr_->getTensorAttrPtr("hencky")->getAttr();
		tensile_limit_ = ensemble_ptr_->getScalarAttrPtr("tensile_limit")->getAttr();
		klog_ = false;
		class_name_ = "StateKLModel";
		//! material parameters for consolidation
		/*
		p_0_ = 16.937;
		Lambda_ = -5.226;
		alpha_ = 9.203;
		*/
		p_0_ = 7.3033;
		alpha_ = 9.4455;
		Lambda_ = 4.7522;
		//! material parameters for critical state
		p_1_ = 12.337;
		p_n_ = -8.5;
		q_n_ = 7.402;
		B_ = 1.168;
		beta_ = 0.5;
		p_i_ = 4.99;

		//! tensile failure pressure
		p_t_ = -0.09;
		tensile_0_ = -p_t_;
		weibull_c_ = 5.23;
		std::default_random_engine generator;
		std::weibull_distribution<double> distribution(weibull_c_, tensile_0_);
		for (int i = 0; i < np_; i++) {
			tensile_limit_[i] = max(0.01, distribution(generator));
		}

		//! parameters for damage
		C1_ = 0.1;
		C2_ = 0.0;
		D1_ = 0.005;
		D2_ = 0.0;
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
		std::pair<double, double> theta = actualTheta(dgamma, p, theta_e, theta_p, damage);
		double theta_ne = theta.first, theta_np = theta.second;
		double damage_n = actualDamage(sigma, dgamma, damage);
		std::pair<Mat3d, Mat3d> sig_h = actualSigmaHencky(h_tr, J, dgamma, theta_ne, theta_np, p, damage_n);
		return yield(sig_h.first, theta_np, damage_n);
	}

	double StateKLModel::pc(const double & theta_p) const {
		//return p_0_ + Lambda_ / alpha_ * (1 - exp(-alpha_ * theta_p));
		return p_0_ + alpha_ * (exp(Lambda_ * theta_p) - 1);
	}

	double StateKLModel::qc(const double & pm) const {
		/*
		if (pm <= p_i_)
			return (p_1_ - pm) / (p_1_ - p_n_) * q_n_;
		else
			return B_ * pow(pm, beta_);
			*/
		return c_[0] + c_[1] * pm + c_[2] * pow(pm, 2) + c_[3] * pow(pm, 3)
			+ c_[4] * pow(pm, 4) + c_[5] * pow(pm, 5);
	}

	double StateKLModel::qf(const double & pc) const {
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
		double b2 = b_[2][0] + b_[2][1] * theta_p + b_[2][2] * pow(theta_p, 2.0) + b_[2][3] * pow(theta_p, 3.0);
		return b0 * theta_e + b1 * pow(theta_e, 2.0) + b2 * pow(theta_e, 3.0);
	}

	double StateKLModel::ddamage(const double& damage, const Mat3d& sigma) const {
		if (damage == 1.0)
			return 0.0;
		double p = -sigma.trace() / 3.0;
		double epf = D1_ * pow(p - p_t_, D2_);
		return sqtwothirds * dev(sigma).norm() / epf;
	}

	//! Function
	double StateKLModel::zeta(const double& p, const double& theta_p, const double& damage) const {
		double p_c = pc(theta_p);
		double q_f = qf(p_c);
		double q_c = qc(0.5 * p_c + 0.5 * p_t_);
		return 4.0 * pow((1 - damage) * q_c + damage * q_f, 2.0) *
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
		double pn = d / J * (dshear(theta_e, theta_p) * pow(J / div, 2.0) * dev(h_tr).squaredNorm() + Jp(theta_e, theta_p));
		Mat3d hencky = h_tr - dgamma * sigma + 1.0 / 3.0 * dgamma * zeta(p, theta_p, damage) * Mat3d::Identity();
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
		Mat3d tp = fb * b_e * fb.transpose();
		//if (tp.determinant() == 0.0) throwException("trialHencky", "determinant is zero!");
		//Mat3d res = tp.log();
		return 0.5 * tp.log().real();// .log();
	}

	Mat3d StateKLModel::trialSigma(const Mat3d& hencky, const double& J, const double& theta_e,
		const double& theta_p, const double& p, const double& damage) const {
		double a = (1.0 - heaviside(-p) * damage) / J;
		Mat3d dev_h = dev(hencky);
		Mat3d dev = a * (2.0 * shear(theta_e, theta_p) * dev_h);
		double p_trial = a * (dshear(theta_e, theta_p) * dev_h.squaredNorm() + Jp(theta_e, theta_p));
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
		double te_tr = -h_tr.trace();
		Mat3d sig_tr = trialSigma(h_tr, J, te_tr, theta_p, p, damage);
		//! check for yield
		if (yield(sig_tr, theta_p, damage) <= epsilon) {
			return Vars(0.0, te_tr, theta_p, sig_tr, h_tr);
		}
		else {
			double dg1 = 0.0;
			double dg2 = 1.0* pow(10.0, -5);
			long cntn = 0;
			double k, min_k = checkYield(h_tr, sigma, J, theta_e, theta_p, p, damage, dg2), min_dg = 0;
			while (true) {
				k = checkYield(h_tr, sigma, J, theta_e, theta_p, p, damage, dg2);
				if (k < 0.0)
					break;
				else {
					if (k < min_k - epsilon) {
						min_k = k;
						min_dg = dg2;
					}
					if (dg2 < 0.3)
						dg2 += pow(10.0, -5);
					else if (dg2 < 0.75)
						dg2 += 5.0 * pow(10.0, -5);
					else if (dg2 < 1.5)
						dg2 += 2.0 * pow(10.0, -4);
					else if (dg2 < 10)
						dg2 += 1.0 * pow(10, -3);
					else dg2 += 1.0 * pow(10, -2);
					if (++cntn == 100000) {
						/*std::string info = "not converging from 0.0 to " + std::to_string((double)(cntn - 3022500) * pow(10.0, -3) + 10.0) 
							+ ", dg2 is: " + std::to_string(min_dg) + "\n yield is: "
							+ std::to_string(min_k);
						throwException("Yield", info);*/
						dg2 = -min_dg;
						break;
					}
				}
			}
			double dg = bisectionSearch(dg1, dg2, h_tr, sigma, J, theta_e, theta_p, p, damage);
			/*
			double dg = 0.0;
			double der = 0.0;
			double tg = 0.0;
			int cntn = 0;
			double cy1, cy2;
			while (true) {
				cy1 = checkYield(h_tr, sigma, J, theta_e, theta_p, p, damage, dg);
				cy2 = checkYield(h_tr, sigma, J, theta_e, theta_p, p, damage, dg + 0.00001);
				if (cy1 != cy2)
					der = (cy2 - cy1) / 0.00001;
				else {
					cy2 = checkYield(h_tr, sigma, J, theta_e, theta_p, p, damage, dg + 0.00005);
					der = (cy2 - cy1) / 0.00005;
				}
				if (abs(cy1 / der) < pow(10, -9) || abs(cy1 / der) > pow(10, 23)) {
					break;
				}
				else {
					tg = dg - cy1 / der;
					if (tg < 0) dg = dg + 0.1;
					else dg = tg;
				}
				if (++cntn == 300) {
					std::string info = "dead loop, d_gamma is: " + std::to_string(dg) + "\n yield is: " + std::to_string(cy1) + 
						"\n der is: " + std::to_string(der) + "\n delta is: " + std::to_string(cy1 / der);
					throwException("yield", info);
				}
			}
			*/
			//std::cout << "delta gamma: " << dg << std::endl;
			double d = actualDamage(sigma, dg, damage);
			std::pair<double, double> theta_a = actualTheta(dg, p, theta_e, theta_p, damage);
			std::pair<Mat3d, Mat3d> sig_h = actualSigmaHencky(h_tr, J, dg, theta_a.first, theta_a.second, p, d);
			return Vars(d, theta_a.first, theta_a.second, sig_h.first, sig_h.second);
		}
	}

	double StateKLModel::bisectionSearch(const double & left, const double & right, const Mat3d & h_tr, 
		const Mat3d & sigma, const double & J, const double & theta_e, 
		const double & theta_p, const double & p, const double & damage) const {
		if (right < 0.0) return -right;
		if (right - left < epsilon) return right;
		double mid = left + (right - left) / 2.0;
		double k = checkYield(h_tr, sigma, J, theta_e, theta_p, p, damage, mid);
		if (k < -epsilon)
			return bisectionSearch(left, mid, h_tr, sigma, J, theta_e, theta_p, p, damage);
		else if (k > epsilon)
			return bisectionSearch(mid, right, h_tr, sigma, J, theta_e, theta_p, p, damage);
		else {
			return right;
		}
	}

}