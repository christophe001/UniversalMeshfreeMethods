
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>
#include "u_timer.h"
//#include "m_neighborhood_test.h"
//#include "m_meshfree_solver.h"
#include "m_collision_plate.h"
using namespace std;
using namespace msl;

int main() {
	/*try {
		NeighborhoodTest nbt(false, 4.0 * Vec3d::Ones(), 0.1);
		nbt.setHorizon(0.35);
		nbt.setSaveParams("F:\\newM4", "TestNeighbor2");
		nbt.runPeri();
		nbt.saveVtk();
	}
	catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
	}*/
	try {
		/*
		DomainConfig cfg(-2.0*Vec3d::Ones(), 2.0*Vec3d::Ones(), 0.2);
		std::shared_ptr<MeshfreeSolver> solver = std::make_shared<MeshfreeSolver>();
		solver->setParams(0.2, 20, false);
		solver->createEnsemble("rectangle", cfg, solver->getJH2AttributeConfig(70.84, 90.16), 
			Vec3d::Ones(), Vec3d::Zero(), Vec3d::UnitZ(), 0, 0.05);
		solver->configRunParams("F:\\newM4", "jh2_model", 0.0001, 2.0);
		solver->configSolver("pd_jh2", "");
		solver->setJH2Params(2);
		solver->run();*/
		
		/*
		Vec3d cmin = Vec3d::Zero();
		Vec3d cmax = Vec3d::Zero();
		for (int i = 0; i < 100; i++) {
			Vec3d temp = i * Vec3d::Ones();
			cmin = cmin.cwiseMin(temp);
			cmax = cmax.cwiseMax(temp);
		}
		std::cout << "cmin: " << cmin.format(PureFmt) << std::endl;
		std::cout << "cmax: " << cmax.format(PureFmt) << std::endl;
		*/
		
		/*
		msl::ObjInfo pt("sphere", Vec3d{ 0.00075, 0.00075, 0.00075 }, "rigid", 0.00008, 
			8500, Vec3d{ 0, 0, -150 }, Vec3d{ 0, 0, 0.00125 });
		msl::ObjInfo targ("rectangle", Vec3d{ 0.02,0.02,0.001 }, "kl", 0.00016, 2200, Vec3d{ 0, 0, 0 });
		DomainConfig cfg(Vec3d{-0.02, -0.02, -0.01}, Vec3d{ 0.02, 0.02, 0.01 }, 0.001);
		std::string folder = "F:\\newM4";
		std::string file = "contact_small";
		*/
		msl::ObjInfo pt("sphere", Vec3d{ 0.0005, 0.0005, 0.0005 }, "rigid", 0.000067,
			7900, Vec3d{ 0, 0, -6900 }, Vec3d{ 0, 0, 0.0031 });
		msl::ObjInfo targ("rectangle", Vec3d{ 0.04,0.04,0.005 }, "kl", 0.00031, 2200, Vec3d{ 0, 0, 0 });
		DomainConfig cfg(Vec3d{ -0.04, -0.04, -0.04 }, Vec3d{ 0.04, 0.04, 0.1 }, 0.0009);
		std::string folder = "F:\\newM4";
		std::string file = "contact_japan_90m_wb_final";

		CollisionPlate cp(pt, targ, cfg, folder, file, 2.5 * pow(10, -10), 10 * pow(10, -6), 50);
		cp.run();
		
	}
	catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
	}
	cin.get();
}