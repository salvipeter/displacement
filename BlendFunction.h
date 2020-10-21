// -*- mode: c++ -*-

#pragma once

#include <geometry.hh>

Geometry::DoubleVector computeBlendFunction(size_t n, const Geometry::TriMesh &mesh,
                                            const Geometry::Point2DVector &parameters,
                                            size_t foot_index, double radius);
