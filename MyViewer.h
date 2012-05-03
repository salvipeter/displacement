#pragma once

#include <string>

#include <QGLViewer/qglviewer.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

class MyViewer : public QGLViewer
{
  Q_OBJECT

public:
  MyViewer(QWidget *parent);
  virtual ~MyViewer();

  inline double getCutoffRatio() const;
  inline void setCutoffRatio(double ratio);
  inline double getMeanMin() const;
  inline void setMeanMin(double min);
  inline double getMeanMax() const;
  inline void setMeanMax(double max);
  bool openMesh(std::string const &filename);

protected:
  virtual void init();
  virtual void draw();

private:
  struct MyTraits : public OpenMesh::DefaultTraits {
    VertexTraits {
      double area;              // total area of the surrounding triangles
      double mean;              // approximated mean curvature
    };
    FaceTraits {
      double area;
    };
  };
  typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;
  typedef MyMesh::Point Vector;

  inline Vector halfedgeVector(MyMesh::HalfedgeHandle const &h) const;
  void updateMeanMinMax();
  void updateMeanCurvature();
  void meanMapColor(double d, double *color) const;

  MyMesh mesh;
  double mean_min, mean_max;
  double cutoff_ratio;
};

#include "MyViewer.hpp"