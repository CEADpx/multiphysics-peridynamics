/*
 * -------------------------------------------
 * Copyright (c) 2021 - 2025 Prashant K. Jha
 * -------------------------------------------
 * https://github.com/CEADpx/multiphysics-peridynamics
 *
 * Distributed under the Boost Software License, Version 1.0. (See accompanying
 * file LICENSE)
 */
#pragma once

#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <libmesh/point.h>
#include "io.h"

namespace geom {


// Base geometry class
class GeomObject {
public:
  std::string d_name;
  std::string d_description;
  std::vector<std::string> d_tags;

public:
  GeomObject(std::string name = "", std::string description = "")
    : d_name(name), d_description(description) {}

  virtual ~GeomObject() = default;

  // Virtual methods
  virtual double volume() const { return 0.; }
     virtual libMesh::Point center() const { return libMesh::Point(); }
   virtual std::pair<libMesh::Point, libMesh::Point> box() const {
     return {libMesh::Point(), libMesh::Point()};
   }
   virtual double inscribedRadius() const { return 0.; }
   virtual double boundingRadius() const { return 0.; }
   
   // Basic geometric queries
   virtual bool isInside(const libMesh::Point& x) const { return false; }
   virtual bool isOutside(const libMesh::Point& x) const { return !isInside(x); }
   virtual bool isNear(const libMesh::Point& x, const double& tol) const { return false; }

  // Getters
  const std::string& getName() const { return d_name; }
  const std::string& getDescription() const { return d_description; }

  // Print helpers
  virtual std::string printStr(int nt = 0, int lvl = 0) const {
    std::ostringstream oss;
    auto tabS = util::io::getTabS(nt);
    oss << tabS << "------- geom::GeomObject --------" << std::endl
        << std::endl;
    oss << tabS << "name = " << d_name << std::endl;
    oss << tabS << "description = " << d_description << std::endl;
    oss << tabS << "tags = ";
    for (auto tag : d_tags) oss << tag << " ";
    oss << std::endl;
    oss << tabS << std::endl;
    return oss.str();
  }

  virtual void print(int nt = 0, int lvl = 0) const { std::cout << printStr(nt, lvl); }
};

// Rectangle implementation
class Rectangle : public GeomObject {
public:
  double d_Lx, d_Ly;
  double d_r;
  libMesh::Point d_x;

public:
  Rectangle(double Lx, double Ly, libMesh::Point x = libMesh::Point(),
           std::string description = "")
    : GeomObject("rectangle", description),
      d_Lx(Lx), d_Ly(Ly),
      d_r(0.5 * std::sqrt(Lx * Lx + Ly * Ly)),
      d_x(x) {}

  double volume() const override { return d_Lx * d_Ly; }
  libMesh::Point center() const override { return d_x; }
  
  std::pair<libMesh::Point, libMesh::Point> box() const override {
    return {
      d_x + libMesh::Point(-0.5 * d_Lx, -0.5 * d_Ly, 0.),
      d_x + libMesh::Point(0.5 * d_Lx, 0.5 * d_Ly, 0.)
    };
  }

  bool isInside(const libMesh::Point& x) const override {
    auto dx = x - d_x;
    return std::abs(dx(0)) <= 0.5 * d_Lx && 
           std::abs(dx(1)) <= 0.5 * d_Ly;
  }

  // Print helpers
  std::string printStr(int nt = 0, int lvl = 0) const override {
    std::ostringstream oss;
    auto tabS = util::io::getTabS(nt);
    oss << tabS << "------- geom::Rectangle --------" << std::endl
        << std::endl;
    oss << tabS << "center = (" << d_x(0) << ", " << d_x(1) << ", " << d_x(2) << ")" << std::endl;
    oss << tabS << "Lx = " << d_Lx << ", Ly = " << d_Ly << std::endl;
    oss << tabS << "Inscribed radius = " << d_r << std::endl;
    oss << tabS << std::endl;
    return oss.str();
  }

  void print(int nt = 0, int lvl = 0) const override { std::cout << printStr(nt, lvl); }
};

// Circle implementation
class Circle : public GeomObject {
public:
  libMesh::Point d_x;
  double d_r;

public:
  Circle(double r, libMesh::Point x = libMesh::Point(),
         std::string description = "")
    : GeomObject("circle", description),
      d_x(x), d_r(r) {}

  double volume() const override { return M_PI * d_r * d_r; }
  libMesh::Point center() const override { return d_x; }
  
  std::pair<libMesh::Point, libMesh::Point> box() const override {
    return {
      d_x + libMesh::Point(-d_r, -d_r, 0.),
      d_x + libMesh::Point(d_r, d_r, 0.)
    };
  }

  bool isInside(const libMesh::Point& x) const override {
    auto dx = x - d_x;
    return dx.norm() <= d_r;
  }

  // Print helpers
  std::string printStr(int nt = 0, int lvl = 0) const override {
    std::ostringstream oss;
    auto tabS = util::io::getTabS(nt);
    oss << tabS << "------- geom::Circle --------" << std::endl
        << std::endl;
    oss << tabS << "center = (" << d_x(0) << ", " << d_x(1) << ", " << d_x(2) << ")" << std::endl;
    oss << tabS << "Radius = " << d_r << std::endl;
    oss << tabS << std::endl;
    return oss.str();
  }

  void print(int nt = 0, int lvl = 0) const override { std::cout << printStr(nt, lvl); }
};

// Cuboid implementation
class Cuboid : public GeomObject {
public:
  double d_Lx, d_Ly, d_Lz;
  double d_r;
  libMesh::Point d_x;

public:
  Cuboid(double Lx, double Ly, double Lz,
         libMesh::Point x = libMesh::Point(),
         std::string description = "")
    : GeomObject("cuboid", description),
      d_Lx(Lx), d_Ly(Ly), d_Lz(Lz),
      d_r(0.5 * std::sqrt(Lx * Lx + Ly * Ly + Lz * Lz)),
      d_x(x) {}

  double volume() const override { return d_Lx * d_Ly * d_Lz; }
  libMesh::Point center() const override { return d_x; }
  
  std::pair<libMesh::Point, libMesh::Point> box() const override {
    return {
      d_x + libMesh::Point(-0.5 * d_Lx, -0.5 * d_Ly, -0.5 * d_Lz),
      d_x + libMesh::Point(0.5 * d_Lx, 0.5 * d_Ly, 0.5 * d_Lz)
    };
  }

  bool isInside(const libMesh::Point& x) const override {
    auto dx = x - d_x;
    return std::abs(dx(0)) <= 0.5 * d_Lx && 
           std::abs(dx(1)) <= 0.5 * d_Ly &&
           std::abs(dx(2)) <= 0.5 * d_Lz;
  }

  // Print helpers
  std::string printStr(int nt = 0, int lvl = 0) const override {
    std::ostringstream oss;
    auto tabS = util::io::getTabS(nt);
    oss << tabS << "------- geom::Cuboid --------" << std::endl
        << std::endl;
    oss << tabS << "center = (" << d_x(0) << ", " << d_x(1) << ", " << d_x(2) << ")" << std::endl;
    oss << tabS << "Lx = " << d_Lx << ", Ly = " << d_Ly << ", Lz = " << d_Lz << std::endl;
    oss << tabS << "Inscribed radius = " << d_r << std::endl;
    oss << tabS << std::endl;
    return oss.str();
  }

  void print(int nt = 0, int lvl = 0) const override { std::cout << printStr(nt, lvl); }
};

// Sphere implementation
class Sphere : public GeomObject {
public:
  libMesh::Point d_x;
  double d_r;

public:
  Sphere(double r, libMesh::Point x = libMesh::Point(),
         std::string description = "")
    : GeomObject("sphere", description),
      d_x(x),
      d_r(r) {}

  double volume() const override { return (4.0/3.0) * M_PI * d_r * d_r * d_r; }
  libMesh::Point center() const override { return d_x; }
  
  std::pair<libMesh::Point, libMesh::Point> box() const override {
    return {
      d_x + libMesh::Point(-d_r, -d_r, -d_r),
      d_x + libMesh::Point(d_r, d_r, d_r)
    };
  }

  bool isInside(const libMesh::Point& x) const override {
    auto dx = x - d_x;
    return dx.norm() <= d_r;
  }

  double inscribedRadius() const override { return d_r; }
  double boundingRadius() const override { return d_r; }

  // Print helpers
  std::string printStr(int nt = 0, int lvl = 0) const override {
    std::ostringstream oss;
    auto tabS = util::io::getTabS(nt);
    oss << tabS << "------- geom::Sphere --------" << std::endl
        << std::endl;
    oss << tabS << "center = (" << d_x(0) << ", " << d_x(1) << ", " << d_x(2) << ")" << std::endl;
    oss << tabS << "Radius = " << d_r << std::endl;
    oss << tabS << std::endl;
    return oss.str();
  }

  void print(int nt = 0, int lvl = 0) const override { std::cout << printStr(nt, lvl); }
};

// GeomData class to store geometry information
class GeomData {
public:
  std::string d_type;
  std::vector<double> d_params;
  std::vector<std::string> d_flags;
  std::shared_ptr<GeomObject> d_geomObject;

public:
  GeomData() = default;

  void setType(const std::string& type) { d_type = type; }
  void setParams(const std::vector<double>& params) { d_params = params; }
  void setFlags(const std::vector<std::string>& flags) { d_flags = flags; }
  void setGeomObject(std::shared_ptr<GeomObject> obj) { d_geomObject = obj; }

  const std::string& getType() const { return d_type; }
  const std::vector<double>& getParams() const { return d_params; }
  const std::vector<std::string>& getFlags() const { return d_flags; }
  std::shared_ptr<GeomObject> getGeomObject() const { return d_geomObject; }

  // Factory method to create geometry objects
     static std::shared_ptr<GeomObject> createGeomObject(
     const std::string& type,
     const std::vector<double>& params,
     const libMesh::Point& center = libMesh::Point()) {
      
    if (type == "rectangle" && params.size() >= 2) {
      return std::make_shared<Rectangle>(params[0], params[1], center);
    }
    else if (type == "circle" && params.size() >= 1) {
      return std::make_shared<Circle>(params[0], center);
    }
    else if (type == "cuboid" && params.size() >= 3) {
      return std::make_shared<Cuboid>(params[0], params[1], params[2], center);
    }
         else if (type == "sphere" && params.size() >= 1) {
       return std::make_shared<Sphere>(params[0], center);
    }
    return nullptr;
  }
};

} // namespace geom 