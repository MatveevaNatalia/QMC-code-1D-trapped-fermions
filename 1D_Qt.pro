TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/Algorithm.cpp \
    src/dens.cpp \
    src/Energy.cpp \
    src/g2D.cpp \
    src/OBDM.cpp \
    src/QMC.cpp \
    src/Statistics.cpp \
    src/Wave_fun.cpp \
    src/main.cpp

HEADERS += \
    src/Algorithm.h \
    src/dens.h \
    src/Energy.h \
    src/g2D.h \
    src/OBDM.h \
    src/Statistics.h \
    src/Wave_fun.h \
    src/qmc.h
