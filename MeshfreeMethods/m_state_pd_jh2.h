/*! \file m_state_pd_jh2.h */

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

#ifndef _M4_STATE_PD_JH2_
#define _M4_STATE_PD_JH2_

#include "m_state_peridynamics.h"

namespace msl {
	class StateBasedPDJH2 : public StateBasedPD {
	protected:
		double K1_, K2_, K3_;
		double A_, B_, C_, M_, N_;
		double T_;
		double D1_, D2_;
		double sigma_hel_, p_hel_;

		double* damage_;
		double* pressure_;
		double* delta_p_;
		//double*	sigma_;
		
	public:
		StateBasedPDJH2(std::shared_ptr<SortEnsemble> sorted, 
			std::shared_ptr<NeighborhoodData> nbh, std::shared_ptr<PeriNeighborData> pbh);

		StateBasedPDJH2(std::shared_ptr<ComputeNeighbor> cpn);

		~StateBasedPDJH2() {}
		
		void computeTau() override;

		void configModelParams(const std::vector<std::vector<double>>& params) override;

		std::string format() const override;
		
		void configParamsK(double k1, double k2 = 0.0, double k3 = 0.0) {
			K1_ = k1; K2_ = k2; K3_ = k3;
		}
		
		void configParamsD(double d1, double d2) {
			D1_ = d1; D2_ = d2;
		}
		
		void configParamsHel(double sigma, double pressure) {
			sigma_hel_ = sigma;
			p_hel_ = pressure;
		}
		
		void configParamsA(double a, double b, double c, double m, double n) {
			A_ = a; B_ = b; C_ = c; M_ = m; N_ = n;
		}


	};
}

#endif // !_M4_STATE_PD_JH2_
