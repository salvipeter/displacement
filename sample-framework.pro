# -*- mode: Makefile -*-

TARGET = sample-framework
CONFIG += c++17 qt opengl debug
QT += gui widgets opengl xml
QMAKE_CXX = clang++

HEADERS = MyWindow.h MyViewer.h MyViewer.hpp BlendFunction.h
SOURCES = MyWindow.cpp MyViewer.cpp main.cpp BlendFunction.cpp

TRANSFINITE = /home/salvi/project/transfinite
INCLUDEPATH += /usr/include/eigen3 $${TRANSFINITE}/src/geom $${TRANSFINITE}/src/transfinite
LIBS *= -lQGLViewer-qt5 -L/usr/lib/OpenMesh -lOpenMeshCore -lGL -lGLU
LIBS += -L$${TRANSFINITE}/release/geom -L$${TRANSFINITE}/release/transfinite -ltransfinite -lgeom

# Optional
# DEFINES += BETTER_MEAN_CURVATURE USE_JET_FITTING
# LIBS += -lCGAL

RESOURCES = sample-framework.qrc
