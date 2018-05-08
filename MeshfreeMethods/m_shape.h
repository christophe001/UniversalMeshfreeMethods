/*! \file m_shape.h */

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

#ifndef _M4_SHAPE_
#define _M4_SHAPE_

#include "m_types.h"
#include "m_obj.h"
#include <vector>
#include <memory>

namespace msl {
	class Shape : public MsObj {
	protected:
		std::string	name_;
		Vec3d		orientation_;
		double		theta_;
		Mat3d		rotation_;
		bool		shape2d_;
	public:
		Shape();
		virtual ~Shape() {}
		void setOrientation(const Vec3d& dir, const double& theta = 0);
		virtual double getFirstDim() const = 0;
		virtual double getVolume() const { return 0.0; }
		virtual void setDims(double a, double b, double c = 0) = 0;
		virtual std::vector<double> getDims() const = 0;
		virtual bool isWithin(const Vec3d& pos) const = 0;
		Vec3d getOrientation() const { return orientation_; }
		double getTheta() const { return theta_; }
		bool shape2D() const { return shape2d_; }
		Mat3d getRotation() const { return rotation_; }
	};

	class Rectangle2D : public Shape {
		double a_, b_;
	public:
		Rectangle2D() { shape2d_ = true; }
		Rectangle2D(double a, double b);
		double getFirstDim() const override { return a_; }
		double getVolume() const override;
		void setDims(double a, double b, double c = 0) override;
		std::vector<double> getDims() const override;
		bool isWithin(const Vec3d& pos) const override;		
	};

	class Circle : public Shape {
		double r_;
	public:
		Circle() { shape2d_ = true; }
		Circle(double r);
		double getFirstDim() const override { return r_; }
		double getVolume() const override;
		void setDims(double a, double b = 0, double c = 0) override;
		std::vector<double> getDims() const override;
		bool isWithin(const Vec3d& pos) const override;
	};

	class Triangle : public Shape {
		double a_, b_;
	public:
		Triangle() { shape2d_ = true; }
		Triangle(double a, double b);
		double getFirstDim() const override { return a_; }
		double getVolume() const override;
		void setDims(double a, double b, double c = 0) override;
		std::vector<double> getDims() const override;
		bool isWithin(const Vec3d& pos) const override;
	};

	class Rectangle : public Shape {
		double a_, b_, c_;
	public:
		Rectangle() { shape2d_ = false; }
		Rectangle(double a, double b, double c);
		double getFirstDim() const override { return a_; }
		double getVolume() const override;
		void setDims(double a, double b, double c ) override;
		std::vector<double> getDims() const override;
		bool isWithin(const Vec3d& pos) const override;
	};

	class Sphere : public Shape {
		double r_;
	public:
		Sphere() { shape2d_ = false; }
		Sphere(double r);
		double getFirstDim() const override { return r_; }
		double getVolume() const  override;
		void setDims(double a, double b = 0, double c = 0) override;
		std::vector<double> getDims() const override;
		bool isWithin(const Vec3d& pos) const override;
	};

	class Cone : public Shape {
		double r_, h_;
	public:
		Cone() { shape2d_ = false; }
		Cone(double r, double h);
		double getFirstDim() const override { return r_; }
		double getVolume() const override;
		void setDims(double a, double b, double c = 0) override;
		std::vector<double> getDims() const override;
		bool isWithin(const Vec3d& pos) const override;
	};

	class Cylinder : public Shape {
		double r_, h_;
	public:
		Cylinder() { shape2d_ = false; }
		Cylinder(double r, double h);
		double getFirstDim() const override { return r_; }
		double getVolume() const override;
		void setDims(double a, double b, double c = 0) override;
		std::vector<double> getDims() const override;
		bool isWithin(const Vec3d& pos) const override;
	};

	class ShapeFactory {
	public:
		static std::shared_ptr<Shape> ShapeBuilder(const std::string& s) {
			if (s == "rectangle2d")
				return std::make_shared<Rectangle2D>();
			else if (s == "cicle")
				return std::make_shared<Circle>();
			else if (s == "triangle")
				return std::make_shared<Triangle>();
			else if (s == "rectangle")
				return std::make_shared<Rectangle>();
			else if (s == "sphere")
				return std::make_shared<Sphere>();
			else if (s == "cone")
				return std::make_shared<Cone>();
			else if (s == "cylinder")
				return std::make_shared<Cylinder>();
			else
				throw std::invalid_argument("unknown shape specifier");
		}
	};
}

#endif // !_ENSEMBLE_CREATOR_