##Templateforautomakeforlibgromos

noinst_LTLIBRARIES = libgmath.la

gincludedir = $(includedir)/gromos++

ginclude_HEADERS = Vec.h \
	Matrix.h \
	Distribution.h \
	WDistribution.h \
	Stat.h \
	Stat.cc \
	StatDisk.h \
	StatDisk.cc \
	Expression.h \
	Correlation.h \
	Physics.h \
	Mesh.h

libgmath_la_SOURCES = Vec.cc \
	Matrix.cc \
	Distribution.cc \
	WDistribution.cc \
	Expression.cc \
	Correlation.cc \
	Physics.cc \
	Mesh.cc

check_PROGRAMS = Vec \
	Matrix \
	Distribution \
	Expression \
	Correlation \
	Mesh

Vec_SOURCES = Vec.t.cc
Matrix_SOURCES = Matrix.t.cc
Distribution_SOURCES = Distribution.t.cc
Expression_SOURCES = Expression.t.cc
Correlation_SOURCES = Correlation.t.cc
Mesh_SOURCES = Mesh.t.cc

LDADD = libgmath.la -lgslcblas -lgsl
