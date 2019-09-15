/*! \file m_types.h*/

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

#ifndef _M4_TYPES_ 
#define _M4_TYPES_

#pragma warning (disable:4244)
#include <Eigen\Dense>
#include <Eigen\Core>
#include <vector>

#ifndef _WITH_OMP_ 
#define _WITH_OMP_
#endif // !_WITH_OMP 

#define _CRT_SECURE_NO_DEPRECATE  
#define _CRT_NONSTDC_NO_DEPRECATE


const double pi = 3.14159265358979323846;

typedef std::vector<std::vector<double>> ModelParams;
typedef std::vector<double> ParamsList;

const ParamsList NOPD{0.2, 0.05, 0.0001,3700};

const ModelParams jh2_param_1{ NOPD, ParamsList{0.93,0,0,0,0.6}, 
		ParamsList{0,0}, ParamsList{2.79, 1.46}, ParamsList{130.95, 0, 0} };

const ModelParams jh2_param_2{ NOPD, ParamsList{0.93,0,0,0,0.6},
		ParamsList{0.005, 1}, ParamsList{ 2.79, 1.46 }, ParamsList{ 130.95, 0, 0 } };

const ModelParams jh2_param_3{ NOPD, ParamsList{0.93,0.31,0,0.6,0.6}, 
		ParamsList{0.005, 1}, ParamsList{ 2.79, 1.46 }, ParamsList{ 130.95, 0, 0 } };

const double rad2deg = 180.0 / pi;

const double sqtwothirds = sqrt(2.0 / 3.0);

const double epsilon = pow(10.0, -9.0);

const Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

const Eigen::IOFormat PureFmt(8, 0, " ", " ");

enum class PeriType {
	kPeriO = 0,
	kPeriX = 1,
	kPeriY = 2,
	kPeriZ = 4,
	kPeriXY = 3,
	kPeriXZ = 5,
	kPeriYZ = 6,
	kPeriXYZ = 7,
};

enum class ScalarAttrType {
	kPid = 0,
	kCellId = 1,
	kIniDensity = 2,
	kDensity = 3,
	kIniPressure = 4,
	kPressure = 5,
	kIniVolume = 6,
	kVolume = 7,
	kEnergy = 8,
	kType = 9,
	kOther = 10,
};

enum class VectorAttrType {
	kIniPosition = 0,
	kIniVelocity = 1,
	kIniAcc = 2,
	kOther = 3,
};

enum class TensorAttrType {
	kIniDeformation = 0,
	kDeformation = 1,
	kOther = 2,
};

typedef Eigen::Vector2d	Vec2d;
typedef Eigen::Vector2i Vec2i;

typedef Eigen::Vector3d Vec3d;
typedef Eigen::Vector3i Vec3i;
typedef Eigen::Vector4d Vec4d;
typedef Eigen::Array3d	Arr3d;
typedef Eigen::Array3i  Arr3i;
typedef Eigen::Array4d  Arr4d;

typedef Eigen::Matrix2d Mat2d;
typedef Eigen::Matrix3d Mat3d;
typedef Eigen::Matrix4d Mat4d;

const Vec3i vec3zero = Vec3i::Zero();
const Vec3i vec3one = Vec3i::Ones();
const Arr3i arr3zero = Arr3i::Zero();
const Arr3i arr3one = Arr3i::Ones();

typedef struct {
	Vec3d	dmin;
	Vec3d	dmax;
	double	dp;
	double  density;
	double  horizon;
	double	dt;
	double	T;
	int		sv_step;
} Params;

typedef struct {
	double	E;			// Young's modulus in GPa
	double	G;			// Fracture energy in J/m^3
	double	horizon;	// horizon in m
	double	m;			// horizon / dp
	double	length;		// length of 2D case in m
	double	width;		// width of 2D case in m
	double	notch;		// length of notch in m
	double	density;    // density of material in kg/m^3
	double  pressure;	// traction pressure in MPa
} Peridm2DParams;

#endif // !_M4_TYPES_ 