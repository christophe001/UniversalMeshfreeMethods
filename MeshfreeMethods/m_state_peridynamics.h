/*! \file m_state_peridynamics.h */

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

#ifndef _M4_STATE_PD_
#define _M4_STATE_PD_

#include "m_lagrangian_compute.h"
#include "m_compute_neighborhood.h"
#include "m_ensemble_creator.h"

namespace msl {
	class StateBasedPD : public LagrangianCompute {
	protected:
		Mat3d*		shape_tensor_;
		Mat3d*		deformation_;
		Mat3d*		deformation_last_;

		Mat3d*		deformation_dot_;
		Mat3d*		rotation_;
		Mat3d*		left_stretch_;
		//Mat3d*		tau_;
		Mat3d*		d_;					//! unrotated deformation rate of tensor D
		Mat3d*		left_stretch_dot_;
		Mat3d*		tau_;				//! cauchy stress

		Vec3d*		init_pos_;

		double*		density_;
		//double*		bond_length_;
		double*		lambda_;			//! lame parameters
		double*		mu_;				//! lame parameters
		int*		nbl_;
		double		horizon_;
		int			np_;
		double		i_epsilon_;			//! influence function parameters			
		double		i_p_;				//! w(xi) = 1/(|xi|+i_epsilon)^i_p_
		double		dv_;				//! delta volume;
		double		dp_;				//! inter particle distance
		double		dt_;				//! dt
		double		rho_;				//! initial density;
		bool		shape_calc_;		//! indicate if calculated shape tensor
			

	public:
		StateBasedPD(std::shared_ptr<SortEnsemble> sorted, 
			std::shared_ptr<NeighborhoodData> nbh, std::shared_ptr<PeriNeighborData> pbh_, double dt = 0);

		StateBasedPD(std::shared_ptr<ComputeNeighbor> cpn);
		
		void configModelParams(const std::vector<std::vector<double>>& params) override;

		std::string format() const override;

		void configParams(double horizon, double dp, double dt, double rho = 1000);
		
		double influence(double d) { return pow(d + i_epsilon_, -i_p_); }
		
		double volumeCorrector(double d);
		
		void addShapeTensor(int i, long it);
		
		void addDeformation(int i, long it);
		
		void addDeformationGradient(int i, long it);
		
		void addForceState(int i, long it);

		void addForceStateNoRot(int i, long it);
		
		void computeShapeTensor();
		
		void computeDeformation();
		
		void computeDeformationDot();

		void computeRotation();
		
		void updateLeftStretch();
		
		void updateDensity();
		
		virtual void computeTau() {}
		
		void computeForces() override;
		
		virtual ~StateBasedPD() {}
	};
}

#endif // !_M4_STATE_PD_


