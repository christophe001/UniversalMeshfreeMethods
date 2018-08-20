/*! \file m_state_kl_model.h */

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

#ifndef _M4_STATE_KL_MODEL_
#define _M4_STATE_KL_MODEL_

#include "m_state_peridynamics.h"

namespace msl {
	class StateKLModel : public StateBasedPD {
	protected:
		bool    klog_;		//! indicate if adopt klog stress rate 
		//*********************************************************************
		//						Material parameters
		//*********************************************************************
		double	D1_, D2_;
		double  C1_, C2_;
		double	p_t_;				//! tensile failure pressure
		double  a_[3][4];
		double  b_[2][4];
		double  p_0_, Lambda_, alpha_;

		//*********************************************************************
		//					Critical state parameters
		//*********************************************************************
		double  p_1_, p_n_, q_n_, B_, beta_, p_i_;

		//*********************************************************************
		//			       Additional material properties
		//*********************************************************************
		double* damage_;		//! damage variable
		double* pressure_;		//! pressure
		double* theta_e_;		//! elastic part of log-volumetric strain 
		double* theta_p_;		//! plastic part of log-volumetric strain
		Mat3d*	hencky_;		//! elastic hencky strain

	public:
		struct Vars {
			double damage;
			double theta_e;
			double theta_p;
			Mat3d hencky;
			Mat3d sigma;
			Vars(double d, double te, double tp, 
				Mat3d s, Mat3d h) {
				damage = d;
				theta_e = te;
				theta_p = tp;
				sigma = s;
				hencky = h;
			}
		};
		void computeTau() override;
		void computeForces() override;
		void adoptKlog() { klog_ = true; }
		void setDamageParams(double c1, double c2, double d1, double d2);
		StateKLModel(std::shared_ptr<SortEnsemble> sorted,
			std::shared_ptr<NeighborhoodData> nbh, std::shared_ptr<PeriNeighborData> pbh, double dt = 0);
		virtual ~StateKLModel() {}
		//! Yield condition
		double yield(const Mat3d& sigma, const double& theta_p, const double& damage) const;
		double checkYield(const Mat3d& h_tr, const Mat3d& sigma, const double& J, const double& theta_e,
			const double& theta_p, const double& p, const double& damage, const double& dgamma) const;
		double pc(const double& theta_p) const;
		double qc(const double& pm) const;
		double qf(const double& pc) const;
		//! Elastic response
		double shear(const double& theta_e, const double& theta_p) const;
		double dshear(const double& theta_e, const double& theta_p) const;
		double Jp(const double& theta_e, const double& theta_p) const;
		//! Damage evolution
		double ddamage(const double& damage, const Mat3d& sigma) const;
		double zeta(const double& p, const double& theta_p, const double& damage) const;
		//! Function
		double heaviside(const double& a) const { return a > 0.0 ? 1.0 : 0.0; }
		Mat3d  dev(const Mat3d& m) const;
		double cbar(const double& b_s, const double& b_t) const;
		double ctilde(const double& b_s, const double& b_t) const;
		//! Return updates based on dgamma
		std::pair<double,double> actualTheta(const double& dgamma, const double& p, const double& theta_e_tr,
			const double& theta_p, const double& damage) const;
		std::pair<Mat3d, Mat3d> actualSigmaHencky(const Mat3d& h_tr, const double& J, const double& dgamma, 
			const double& theta_e, const double& theta_p, const double& p, const double& damage) const;
		double actualDamage(const Mat3d& sigma, const double& dgamma, const double& damage) const;
		//! trial functions
		Mat3d trialHencky(const Mat3d& f_last, const Mat3d& f, const Mat3d& b_e) const;
		Mat3d trialSigma(const Mat3d& hencky, const double& J, const double& theta_e,
			const double& theta_p, const double& p, const double& damage) const;
		//! Lagrangian multiplier calculation
		Vars dgamma(const Mat3d& sigma, const Mat3d& hencky, const Mat3d& f_last, const Mat3d& f,
			const double& theta_e, const double& theta_p, const double& damage) const;
	};
}

#endif // !_M4_STATE_KL_MODEL_

