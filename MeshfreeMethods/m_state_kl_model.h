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
			std::shared_ptr<NeighborhoodData> nbh, std::shared_ptr<PeriNeighborData> pbh);
		virtual ~StateKLModel() {}
		//! Yield condition
		double yield(Mat3d sigma, double theta_p, double damage) const;
		double checkYield(Mat3d h_tr, Mat3d sigma, double J, double theta_e,
			double theta_p, double p, double damage, double dgamma) const;
		double pc(double theta_p) const;
		double qc(double pm) const;
		double qf(double pc) const;
		//! Elastic response
		double shear(double theta_e, double theta_p) const;
		double dshear(double theta_e, double theta_p) const;
		double Jp(double theta_e, double theta_p) const;
		//! Damage evolution
		double ddamage(double damage, Mat3d sigma) const;
		double zeta(double p, double theta_p, double damage) const;
		//! Function
		double heaviside(double a) const { return a > 0.0 ? 1.0 : 0.0; }
		Mat3d  dev(Mat3d m) const;
		double cbar(double b_s, double b_t) const;
		double ctilde(double b_s, double b_t) const;
		//! Return updates based on dgamma
		std::pair<double,double> actualTheta( double dgamma, double p, double theta_e_tr,
			double theta_p, double damage) const;
		std::pair<Mat3d, Mat3d> actualSigmaHencky(Mat3d h_tr, double J, double dgamma, double theta_e, 
			double theta_p, double p, double damage) const;
		double actualDamage(Mat3d sigma, double dgamma, double damage) const;
		//! trial functions
		Mat3d trialHencky(Mat3d f_last, Mat3d f, Mat3d b_e) const;
		Mat3d trialSigma(Mat3d hencky, double J, double theta_e, 
			double theta_p, double p, double damage) const;
		//! Lagrangian multiplier calculation
		Vars dgamma(Mat3d sigma, Mat3d hencky, Mat3d f_last, Mat3d f, double J, 
			double theta_e, double theta_p, double damage) const;
	};
}

#endif // !_M4_STATE_KL_MODEL_

