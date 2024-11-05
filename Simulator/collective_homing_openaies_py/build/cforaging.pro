######################
# Evorobot* makefile #
######################

TEMPLATE = app
TARGET = cforaging 

QT += widgets

# Python library

INCLUDEPATH += /usr/include/python2.7
LIBS += -L"/usr/lib/python2.7" -lpython2.7

# GSL library, lib and include path
INCLUDEPATH += /usr/include/gsl
LIBS += -L"/usr/lib/x86_64-linux-gnu" -lgsl -lgslcblas

# Eigen and unsupported include path 
INCLUDEPATH += /home/paolo/aggregation/collective_homing_openaies_py/src

HEADERS += ../src/main.h 
HEADERS += ../src/cforaging.h 
HEADERS += ../src/evonet.h 
HEADERS += ../src/evoopenaies.h 
HEADERS += ../src/robot-env.h 
HEADERS += ../src/utilities.h

SOURCES += ../src/main.cpp 
SOURCES += ../src/cforaging.cpp 
SOURCES += ../src/evonet.cpp 
SOURCES += ../src/evoopenaies.cpp
SOURCES += ../src/robot-env.cpp
SOURCES += ../src/utilities.cpp

