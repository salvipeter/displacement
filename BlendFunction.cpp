#include <Eigen/Core>

#include "BlendFunction.h"

using namespace Geometry;

static double trapezoid(const std::function<double(double)> &f, double a, double b,
                        size_t n, double s) {
  if (n == 1)
    return (f(a) + f(b)) / 2.0 * (b - a);
  double k = std::pow(2, n - 2);
  double h = (b - a) / k;
  double sum = 0;
  for (double x = a + h / 2.0; x <= b; x += h)
    sum += f(x);
  return (s + h * sum) / 2.0;
}

static double simpson(const std::function<double(double)> &f, double a, double b,
                      size_t iterations = 100, double epsilon = 1.0e-7) {
  double s, s_prev = std::numeric_limits<double>::lowest(), st_prev = s_prev;
  for (size_t i = 1; i <= iterations; ++i) {
    double st = trapezoid(f, a, b, i, st_prev);
    s = (4.0 * st - st_prev) / 3.0;
    if (i > 5 && (std::abs(s - s_prev) < epsilon * std::abs(s_prev) || (s == 0 && s_prev == 0)))
      break;
    s_prev = s;
    st_prev = st;
  }
  return s;
}

static double erbsBlend(double t) {
  size_t iterations = 20;
  constexpr double Sd = 1.6571376796460222;
  auto phi = [](double s) { return std::exp(-std::pow(s - 0.5, 2) / (s * (1 - s))); };
  return Sd * simpson(phi, 0, t, iterations);
}

[[maybe_unused]]
static Point2D intersectCircle(double r, const Point2D &a, const Point2D &b) {
  auto d = b - a;
  auto dr2 = d.normSqr();
  auto D = a[0] * b[1] - a[1] * b[0];
  auto R = std::sqrt(r * r * dr2 - D * D);
  auto sign = d[1] < 0 ? -1.0 : 1.0;
  auto x1 = (D * d[1] + sign * d[0] * R) / dr2;
  auto x2 = (D * d[1] - sign * d[0] * R) / dr2;
  auto y1 = (-D * d[0] + std::abs(d[1]) * R) / dr2;
  auto y2 = (-D * d[0] - std::abs(d[1]) * R) / dr2;
  Point2D p1(x1, y1), p2(x2, y2);
  return (p1 - a) * (b - a) < 0 ? p2 : p1;
}

static Point2D intersectLines(const Point2D &p1, const Point2D &p2,
                              const Point2D &q1, const Point2D &q2) {
  auto ap = p1, ad = p2 - p1, bp = q1, bd = q2 - q1;
  double a = ad * ad, b = ad * bd, c = bd * bd;
  double d = ad * (ap - bp), e = bd * (ap - bp);
  if (a * c - b * b < 1.0e-7)
    return ap;
  double s = (b * e - c * d) / (a * c - b * b);
  return ap + ad * s;
}

[[maybe_unused]]
static Point2D intersectPoly(size_t n, const Point2D &a, const Point2D &b) {
  Point2DVector poly;
  for (size_t i = 0; i < n; ++i) {
    double phi = 2 * i * M_PI / n;
    poly.emplace_back(std::cos(phi), std::sin(phi));
  }
  auto dmin = std::numeric_limits<double>::max();
  Point2D pmin;
  for (size_t i = 0; i < n; ++i) {
    const auto &p1 = poly[i];
    const auto &p2 = poly[(i+1)%n];
    auto x = intersectLines(a, b, p1, p2);
    if ((x - a) * (b - a) < 0)
      continue;
    double d = (x - a).norm();
    if (d < dmin) {
      dmin = d;
      pmin = x;
    }
  }
  return pmin;
}

[[maybe_unused]]
static double displacementBlend(size_t n, const Point2D &footpoint, const Point2D &uv, double r) {
  double d = (uv - footpoint).norm();
  if (d < 1e-5)
    return 1.0;
  double dmax = (intersectCircle(std::sin(M_PI / 2 - M_PI / n), footpoint, uv) - footpoint).norm();
  // double dmax = (intersectPoly(n, footpoint, uv) - footpoint).norm();
  if (d - dmax > -1e-5)
    return 0.0;
  double alpha = erbsBlend(d / dmax);
  double h = std::min(d / ((1 - alpha) * r + alpha * dmax), 1.0);
  return 1 - erbsBlend(h);
}

DoubleVector computeBlendFunction(size_t n, const TriMesh &mesh, const Point2DVector &parameters,
                                  size_t foot_index, double radius) {
  DoubleVector result;

  // Old version
  const auto &foot = parameters[foot_index];
  result.reserve(parameters.size());
  for (const auto &uv : parameters)
    result.push_back(displacementBlend(n, foot, uv, radius));

  return result;
}
