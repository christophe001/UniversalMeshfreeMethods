/*! \file m_model_vtk.cpp */

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

#include "m_model_vtk.h"
#include <windows.h>
#include <fstream>
#include <iostream>

namespace msl {
	ModelVtk::ModelVtk() {
		class_name_ = "ModelVtk";
	}
	void ModelVtk::config(std::string proj, Vec3d p, std::string targ, Vec3d t, double dp, double disp) {
		proj_ = proj;
		p_dims_ = p;
		targ_ = targ;
		t_dims_ = t;
		dp_ = dp;
		disp_ = disp;
	}

	void ModelVtk::configIO(std::string folder, std::string file) {
		folder_ = folder;
		file_ = file;
		std::string inner_folder = folder_ + "\\" + file_;
		CreateDirectory(std::wstring(folder_.begin(), folder_.end()).c_str(), NULL);
		CreateDirectory(std::wstring(inner_folder.begin(), inner_folder.end()).c_str(), NULL);
	}

	void ModelVtk::saveVtk() {
		Vec3d p_off{ 0, 0, disp_ };
		auto p_creator = std::make_shared<EnsembleCreator>(proj_, p_off);
		auto t_creator = std::make_shared<EnsembleCreator>(targ_, Vec3d::Zero());
		p_creator->setDims(p_dims_);
		p_creator->setDpActual(dp_);
		p_creator->create();
		int p_np = p_creator->getEnsemble()->getNp();
		std::cout << p_np << std::endl;
		Vec3d* p_pos = p_creator->getEnsemble()->getPos();

		t_creator->setDims(t_dims_);
		t_creator->setDpActual(dp_);
		t_creator->create();
		int t_np = t_creator->getEnsemble()->getNp();
		std::cout << t_np << std::endl;
		Vec3d* t_pos = t_creator->getEnsemble()->getPos();

		std::string inner_folder = folder_ + "\\" + file_;
		std::string vtkname = file_ + ".vtk";
		std::ofstream file(inner_folder + "\\" + vtkname, std::ios::binary);
		int np = p_np + t_np;
		file << "# vtk DataFile Version 3.0\n"
			<< "vtk output\n"
			<< "ASCII\n"
			<< "DATASET POLYDATA\n"
			<< "POINTS " << np << " double\n";
		for (int i = 0; i < p_np; i++) {
			file << p_pos[i].format(PureFmt) << "\n";
		}
		for (int i = 0; i < t_np; i++) {
			file << t_pos[i].format(PureFmt) << "\n";
		}

		file << "VERTICES " << np << " " << np * 2 << "\n";
		for (int i = 0; i < np; i++)
			file << "1 " << i << "\n";
		file << "POINT_DATA " << np << "\n";
		file << "SCALARS Object double\n";
		file << "LOOKUP_TABLE default\n";
		for (int i = 0; i < p_np; i++)
			file << 2 << "\n";
		for (int i = 0; i < t_np; i++)
			file << 1 << "\n";
		file.close();
	}
}