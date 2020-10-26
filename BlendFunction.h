// -*- mode: c++ -*-

#pragma once

#include <geometry.hh>

Geometry::DoubleVector computeBlendFunction(size_t n, const Geometry::TriMesh &mesh,
                                            const Geometry::Point2DVector &parameters,
                                            const std::function<bool(size_t i)> &on_edge,
                                            size_t foot_index, double radius);
