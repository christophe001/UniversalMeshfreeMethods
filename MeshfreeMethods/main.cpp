#include "m_types.h"
#include <iostream>
#include "m_shape.h"
#include <memory>
#include <unsupported/Eigen/MatrixFunctions>
#include "u_timer.h"
#include "m_neighborhood_test.h"
#include "m_meshfree_solver.h"

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
		DomainConfig cfg(-2.0*Vec3d::Ones(), 2.0*Vec3d::Ones(), 0.2);
		std::shared_ptr<MeshfreeSolver> solver = std::make_shared<MeshfreeSolver>();
		solver->setParams(0.2, 20, false);
		solver->createEnsemble("rectangle", cfg, solver->getJH2AttributeConfig(70.84, 90.16), 
			Vec3d::Ones(), Vec3d::Zero(), Vec3d::UnitZ(), 0, 0.05);
		solver->configRunParams("F:\\newM4", "jh2_model", 0.0001, 2.0);
		solver->configSolver("pd_jh2", "");
		solver->setJH2Params(2);
		solver->run();

	}
	catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
	}
	cin.get();
}