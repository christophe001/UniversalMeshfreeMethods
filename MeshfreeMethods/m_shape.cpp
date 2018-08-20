/*! \file m_shape.cpp */

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

#include "m_shape.h"
#include <Eigen/Geometry>

namespace msl {

	//==============================================================================
	/// Shape
	//==============================================================================
	Shape::Shape() : name_("shape") {
		class_name_ = "Shape";
		orientation_ = Vec3d::UnitZ();
		rotation_ = Mat3d::Identity();
		theta_ = 0.0;
		shape2d_ = true;
	}


	//==============================================================================
	/// ShapeBuilder
	/// dir: direction of shape axis, theta: angle rotate aound axis(counter-clockwise)
	//==============================================================================
	void Shape::setOrientation(const Vec3d & dir, const double & theta) {
		if (dir.norm() == 0)
			throwException("setOrientation", "direction norm is zero");
		orientation_ = dir / dir.norm();
		theta_ = theta;
		Mat3d r1, r2;
		if (orientation_ == Vec3d::UnitZ())
			r1 = Mat3d::Identity();
		else {
			Vec3d axis1 = Vec3d::UnitZ().cross(orientation_);
			double ang1 = asin(axis1.norm());
			Vec3d axis = axis1 / axis1.norm();
			r1 = Eigen::AngleAxisd(ang1, axis);
		}
		r2 = Eigen::AngleAxisd(theta_, orientation_);
		rotation_ = r1 * r2;
	}

	//==============================================================================
	/// Rectangle2D
	//==============================================================================
	Rectangle2D::Rectangle2D(double a, double b) {
		if (a <= 0 || b <= 0)
			throwException("Rectangle2D", "Illegal dimensions");
		a_ = a;
		b_ = b;
		shape2d_ = true;
		name_ = "rectangle_2d";
	}

	double Rectangle2D::getVolume() const {
		return a_ * b_;
	}

	void Rectangle2D::setDims(double a, double b, double c) {
		if (a <= 0 || b <= 0)
			throwException("Rectangle2D", "Illegal dimensions");
		a_ = a;
		b_ = b;
	}

	std::vector<double> Rectangle2D::getDims() const {
		std::vector<double> dims(3, 0.0);
		dims[0] = a_;
		dims[1] = 0.0;
		dims[2] = b_;
		return dims;
	}

	bool Rectangle2D::isWithin(const Vec3d& pos) const {
		return abs(pos[0]) < a_ / 2.0 && abs(pos[2]) < b_ / 2.0 && abs(pos[1]) < epsilon;
	}

	//==============================================================================
	/// Circle
	//==============================================================================
	Circle::Circle(double r) {
		if (r <= 0)
			throwException("Circle", "Illegal dimensions");
		r_ = r;
		shape2d_ = true;
		name_ = "circle";
	}

	double Circle::getVolume() const {
		return pi * r_ * r_;
	}

	void Circle::setDims(double a, double b, double c) {
		if (a <= 0)
			throwException("Circle", "Illegal dimensions");
		r_ = a;
	}

	std::vector<double> Circle::getDims() const {
		std::vector<double> dims(3, 0.0);
		dims[0] = 2.0 * r_;
		dims[1] = 0.0;
		dims[2] = 2.0 * r_;
		return dims;
	}
	
	bool Circle::isWithin(const Vec3d& pos) const {
		return pos.norm() < r_ && abs(pos[1]) < epsilon;
	}

	//==============================================================================
	/// Triangle
	//==============================================================================
	Triangle::Triangle(double a, double b) {
		if (a <= 0 || b <= 0)
			throwException("Triangle", "Illegal dimensions");
		a_ = a;
		b_ = b;
		shape2d_ = true;
		name_ = "triangle";
	}

	double Triangle::getVolume() const {
		return a_ * b_ / 2.0;
	}

	void Triangle::setDims(double a, double b, double c) {
		if (a <= 0 || b <= 0)
			throwException("Triangle", "Illegal dimensions");
		a_ = a;
		b_ = b;
	}

	std::vector<double> Triangle::getDims() const {
		std::vector<double> dims(3, 0.0);
		dims[0] = a_; 
		dims[1] = 0.0;
		dims[2] = b_;
		return dims;
	}

	bool Triangle::isWithin(const Vec3d& pos) const {
		return abs(pos[1]) < epsilon && abs(pos[2]) < b_/2.0 &&  
			abs(pos[0]) < a_ / (2.0 * b_) * (pos[2] + b_ / 2.0);

	}

	//==============================================================================
	/// Rectangle
	//==============================================================================
	Rectangle::Rectangle(double a, double b, double c) {
		if (a <= 0 || b <= 0 || c <= 0)
			throwException("Rectangle", "Illegal dimensions");
		a_ = a;
		b_ = b;
		c_ = c;
		shape2d_ = false;
		name_ = "rectangle";
	}

	double Rectangle::getVolume() const {
		return a_ * b_ * c_;
	}

	void Rectangle::setDims(double a, double b, double c) {
		if (a <= 0 || b <= 0 || c <= 0)
			throwException("Rectangle", "Illegal dimensions");
		a_ = a;
		b_ = b;
		c_ = c;
	}

	std::vector<double> Rectangle::getDims() const {
		std::vector<double> dims(3, 0.0);
		dims[0] = a_; 
		dims[1] = b_;
		dims[2] = c_;
		return dims;
	}

	bool Rectangle::isWithin(const Vec3d& pos) const {
		return abs(pos[0]) < a_ / 2.0 && abs(pos[1]) < b_/2.0 && abs(pos[2]) < c_/2.0;
	}

	//==============================================================================
	/// Sphere
	//==============================================================================
	Sphere::Sphere(double r) {
		if (r < 0)
			throwException("Sphere", "Illegal dimensions");
		r_ = r;
		shape2d_ = false;
		name_ = "sphere";
	}

	double Sphere::getVolume() const {
		return 4.0 / 3.0*pi*pow(r_, 3.0);
	}

	void Sphere::setDims(double a, double b, double c) {
		if (a < 0)
			throwException("Sphere", "Illegal dimensions");
		r_ = a;
	}

	std::vector<double> Sphere::getDims() const {
		std::vector<double> dims(3, 2.0 * r_);
		return dims;
	}

	bool Sphere::isWithin(const Vec3d& pos) const {
		return pos.norm() < r_;
	}

	//==============================================================================
	/// Cone
	//==============================================================================
	Cone::Cone(double r, double h) {
		if (r <= 0 || h <= 0)
			throwException("Cone", "Illegal dimensions");
		r_ = r;
		h_ = h;
		shape2d_ = false;
		name_ = "cone";
	}

	double Cone::getVolume() const {
		return pi / 3.0 * r_ * r_ * h_;
	}

	void Cone::setDims(double r, double h, double c) {
		if (r <= 0 || h <= 0)
			throwException("Cone", "Illegal dimensions");
		r_ = r;
		h_ = h;
	}

	std::vector<double> Cone::getDims() const {
		std::vector<double> dims(3, 0.0);
		dims[0] = 2.0 * r_;
		dims[1] = 2.0 * r_;
		dims[2] = h_;
		return dims;
	}

	bool Cone::isWithin(const Vec3d& pos) const {
		return abs(pos[2]) < h_/2.0 && sqrt(pos[0] * pos[0] + pos[1] * pos[1]) < r_ / h_*(pos[2]+ h_/2.0);
	}

	//==============================================================================
	/// Cylinder
	//==============================================================================
	Cylinder::Cylinder(double r, double h) {
		if (r <= 0 || h <= 0)
			throwException("Cylinder", "Illegal dimensions");
		r_ = r;
		h_ = h;
		shape2d_ = false;
		name_ = "cylinder";
	}

	double Cylinder::getVolume() const {
		return pi * r_ * r_ * h_;
	}

	void Cylinder::setDims(double r, double h, double c) {
		if (r <= 0 || h <= 0)
			throwException("Cylinder", "Illegal dimensions");
		r_ = r;
		h_ = h;
	}

	std::vector<double> Cylinder::getDims() const {
		std::vector<double> dims(3, 0.0);
		dims[0] = 2.0 * r_;
		dims[1] = 2.0 * r_;
		dims[2] = h_;
		return dims;
	}

	bool Cylinder::isWithin(const Vec3d& pos) const {
		return abs(pos[2]) < h_ / 2.0 && sqrt(pos[0] * pos[0] + pos[1] * pos[1]) < r_ /2.0;
	}
}