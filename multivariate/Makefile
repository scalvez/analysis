#
# Collie Makefile  V00-04-10
#

# This is the configurable region
# Point these variables to the corrsponding place on your system
# The current setup only requires one to specify the location of
# the files needed for ROOT.

ROOT := ${ROOTSYS}
ROOTINC := ${ROOT}/include/root
ROOTLIB := ${ROOT}/lib/root
ROOTBIN := ${ROOT}/bin

#Pointers to collie areas
IO_INC := ${COLLIESYS}/io/include/
IO_SRC := ${COLLIESYS}/io/src/
LIM_INC := ${COLLIESYS}/limit/include/
LIM_SRC := ${COLLIESYS}/limit/src/
MIN_INC := ${COLLIESYS}/minuit/include/
MIN_SRC := ${COLLIESYS}/minuit/src/
MIN_LIB := ${COLLIESYS}/minuit/lib/
CLHEP_INC  := ${COLLIESYS}
CLHEP_SRC  := ${COLLIESYS}/CLHEP/src
CLHEP_LIB  := ${COLLIESYS}/CLHEP/lib
COLLIELIB := ${COLLIESYS}/lib/
COLLIEEXMP := ${COLLIESYS}/examples/

# end of the configurable region

# Directories to search for header files
SEARCHDIRS :=  -I. -I${ROOTINC} -I${IO_INC} -I${LIM_INC} -I${CLHEP_INC} -I${MIN_INC}


# variables

LINKER       := gcc
DEPENDFLAGS  := -g -O3 -ffast-math -fno-math-errno -Wall ${SEARCHDIRS} -fPIC -Df2cFortran -s
TOUCHHEADERS := ${MYCODEDIR}/*.h


# FORTRAN

F77     := gfortran
FFLAGS  = -fPIC

# C

CC     := gcc
CFLAGS  = ${DEPENDFLAGS}

# C++

CXX      := g++
CXXFLAGS  = ${DEPENDFLAGS}

%.o : %.f
	${F77} ${FFLAGS} -c $< -o $@
%.o : %.c
	${CC} ${CFLAGS} -c $< -o $@
%.o : %.cc
	${CXX} ${CPPFLAGS} -c $< -o $@ ${CXXFLAGS}
%.o : %.cpp
	${CXX} ${CPPFLAGS} -c $< -o $@ ${CXXFLAGS}
%.o : %.cxx
	${CXX} ${CPPFLAGS} -c $< -o $@ ${CXXFLAGS} -I${ROOTINC}/root

# C preprocessor (C, C++, FORTRAN)

CPPFLAGS =

# linker
LOADLIBS := -lgfortran -lm -lstdc++ -ldl
LIMLIBS := ${MIN_LIB}/libCollieMinuit.a -lCollieCLHEP -lnss_nis

ROOTLIBS := -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread
EXAMPLELIBS := ${ROOTLIBS}  -lm -ldl -rdynamic -lCollieIO -lCollieLimit -lCollieCLHEP
LDFLAGS    = -L${ROOTLIB} -L${COLLIELIB} -L${MIN_LIB} -L${CLHEP_LIB}

LIMITSRCS := ${LIM_SRC}/SigBkgdDist.cc ${LIM_SRC}/CLpoint.cc ${LIM_SRC}/CLcompute.cc
LIMITSRCS += ${LIM_SRC}/FitTest.cc ${LIM_SRC}/Loader.cc ${LIM_SRC}/CLfit2.cc ${LIM_SRC}/CLfitSS.cc
LIMITSRCS += ${LIM_SRC}/CLfast.cc ${LIM_SRC}/CLsyst.cc ${LIM_SRC}/CLfit.cc
LIMITSRCS += ${LIM_SRC}/CollieLoader.cc ${LIM_SRC}/FileSetLoader.cc
LIMITSRCS += ${LIM_SRC}/ProfileLH.cc ${LIM_SRC}/ProfileLH_2D.cc ${LIM_SRC}/CrossSectionLimit_FastApprox.cc
LIMITSRCS += ${LIM_SRC}/CrossSectionLimit.cc ${LIM_SRC}/ThreeSigmaEvidence.cc
LIMITSRCS += ${LIM_SRC}/ExclusionLimit.cc ${LIM_SRC}/CrossSectionCalc.cc
LIMITSRCS += ${LIM_SRC}/timeBasedSeed.cc ${LIM_SRC}/CrossSectionCalc2D.cc
LIMITOBJS :=$(patsubst %.cpp,%.o,$(patsubst %.cxx,%.o,$(patsubst %.cc,%.o,${LIMITSRCS})))

LIMITHEADERS := ${LIM_INC}/CLpoint.hh

MINUITSRCS := ${MIN_SRC}/intrac.c ${MIN_SRC}/mncler.f ${MIN_SRC}/mncros.f ${MIN_SRC}/mnemat.f
MINUITSRCS += ${MIN_SRC}/mnfixp.f ${MIN_SRC}/mnhess.f ${MIN_SRC}/mnintr.f ${MIN_SRC}/mnmnos.f
MINUITSRCS += ${MIN_SRC}/mnpint.f ${MIN_SRC}/mnrazz.f ${MIN_SRC}/mnscan.f ${MIN_SRC}/mnstat.f
MINUITSRCS += ${MIN_SRC}/mnvert.f ${MIN_SRC}/minuit.f ${MIN_SRC}/mncntr.f ${MIN_SRC}/mncuve.f
MINUITSRCS += ${MIN_SRC}/mnerrs.f ${MIN_SRC}/mnfree.f ${MIN_SRC}/mnimpr.f ${MIN_SRC}/mnlims.f
MINUITSRCS += ${MIN_SRC}/mnmnot.f ${MIN_SRC}/mnplot.f ${MIN_SRC}/mnread.f ${MIN_SRC}/mnseek.f
MINUITSRCS += ${MIN_SRC}/mnstin.f ${MIN_SRC}/mnwarn.f ${MIN_SRC}/mnamin.f ${MIN_SRC}/mncomd.f
MINUITSRCS += ${MIN_SRC}/mnderi.f ${MIN_SRC}/mneval.f ${MIN_SRC}/mngrad.f ${MIN_SRC}/mninex.f
MINUITSRCS += ${MIN_SRC}/mnline.f ${MIN_SRC}/mnparm.f ${MIN_SRC}/mnpout.f ${MIN_SRC}/mnrn15.f
MINUITSRCS += ${MIN_SRC}/mnset.f ${MIN_SRC}/mntiny.f ${MIN_SRC}/mnwerr.f ${MIN_SRC}/mnbins.f
MINUITSRCS += ${MIN_SRC}/mncont.f ${MIN_SRC}/mndxdi.f ${MIN_SRC}/mnexcm.f ${MIN_SRC}/mnhelp.f
MINUITSRCS += ${MIN_SRC}/mninit.f ${MIN_SRC}/mnmatu.f ${MIN_SRC}/mnpars.f ${MIN_SRC}/mnprin.f
MINUITSRCS += ${MIN_SRC}/mnrset.f ${MIN_SRC}/mnseti.f ${MIN_SRC}/mnunpt.f  ${MIN_SRC}/stand.f
MINUITSRCS += ${MIN_SRC}/mncalf.f ${MIN_SRC}/mncrck.f ${MIN_SRC}/mneig.f ${MIN_SRC}/mnexin.f
MINUITSRCS += ${MIN_SRC}/mnhes1.f ${MIN_SRC}/mninpu.f ${MIN_SRC}/mnmigr.f ${MIN_SRC}/mnpfit.f
MINUITSRCS += ${MIN_SRC}/mnpsdf.f ${MIN_SRC}/mnsave.f ${MIN_SRC}/mnsimp.f ${MIN_SRC}/mnvers.f
MINUITOBJS :=$(patsubst %.f,%.o,$(patsubst %.c,%.o,${MINUITSRCS}))

CLHEPSRCS := ${CLHEP_SRC}/DoubConv.cc ${CLHEP_SRC}/Hurd160Engine.cc ${CLHEP_SRC}/RandChiSquare.cc
CLHEPSRCS += ${CLHEP_SRC}/RandGeneral.cc ${CLHEP_SRC}/RanecuEngine.cc ${CLHEP_SRC}/DRand48Engine.cc
CLHEPSRCS += ${CLHEP_SRC}/Hurd288Engine.cc ${CLHEP_SRC}/RandEngine.cc ${CLHEP_SRC}/RandLandau.cc
CLHEPSRCS += ${CLHEP_SRC}/Ranlux64Engine.cc ${CLHEP_SRC}/DualRand.cc ${CLHEP_SRC}/JamesRandom.cc
CLHEPSRCS += ${CLHEP_SRC}/RandExponential.cc ${CLHEP_SRC}/Random.cc ${CLHEP_SRC}/RanluxEngine.cc
CLHEPSRCS += ${CLHEP_SRC}/EngineFactory.cc ${CLHEP_SRC}/MTwistEngine.cc ${CLHEP_SRC}/RandFlat.cc
CLHEPSRCS += ${CLHEP_SRC}/RandomEngine.cc ${CLHEP_SRC}/RanshiEngine.cc ${CLHEP_SRC}/engineIDulong.cc
CLHEPSRCS += ${CLHEP_SRC}/NonRandomEngine.cc ${CLHEP_SRC}/RandGamma.cc ${CLHEP_SRC}/RandPoisson.cc
CLHEPSRCS += ${CLHEP_SRC}/StaticRandomStates.cc ${CLHEP_SRC}/erfQ.cc ${CLHEP_SRC}/RandBinomial.cc
CLHEPSRCS += ${CLHEP_SRC}/RandGauss.cc ${CLHEP_SRC}/RandPoissonQ.cc ${CLHEP_SRC}/TripleRand.cc
CLHEPSRCS += ${CLHEP_SRC}/flatToGaussian.cc ${CLHEP_SRC}/RandBit.cc ${CLHEP_SRC}/RandGaussQ.cc
CLHEPSRCS += ${CLHEP_SRC}/RandPoissonT.cc ${CLHEP_SRC}/gammln.cc ${CLHEP_SRC}/RandBreitWigner.cc
CLHEPSRCS += ${CLHEP_SRC}/RandGaussT.cc ${CLHEP_SRC}/RandStudentT.cc
CLHEPOBJS :=$(patsubst %.cc,%.o,$(patsubst %.c,%.o,${CLHEPSRCS}))

IOSRCS := ${IO_SRC}/CollieDistribution.cc ${IO_SRC}/CollieEventList.cc
IOSRCS += ${IO_SRC}/CollieChannel.cc ${IO_SRC}/CollieIterator.cc
IOSRCS += ${IO_SRC}/CollieHistogram.cc ${IO_SRC}/CollieQuickKEYS.cc
IOSRCS += ${IO_SRC}/CollieIOFile.cc ${IO_SRC}/CollieXsec.cc
IOSRCS += ${IO_SRC}/CollieMasspoint.cc ${IO_SRC}/CDF_DZero_Distribution.cc
IOSRCS += ${IO_SRC}/CDF_DZero_IOpoint.cc ${IO_SRC}/CDF_DZero_IOfile.cc
IOSRCS += ${IO_SRC}/collieio_dict.cxx
IOOBJS :=$(patsubst %.cpp,%.o,$(patsubst %.cxx,%.o,$(patsubst %.cc,%.o,${IOSRCS})))

IOHEADERS := ${IO_INC}/CollieDistribution.hh ${IO_INC}/CollieEventList.hh
IOHEADERS += ${IO_INC}/CollieChannel.hh ${IO_INC}/CollieIterator.hh
IOHEADERS += ${IO_INC}/CollieHistogram.hh ${IO_INC}/CollieQuickKEYS.hh
IOHEADERS += ${IO_INC}/CollieIOFile.hh ${IO_INC}/CollieXsec.hh
IOHEADERS += ${IO_INC}/CollieMasspoint.hh ${IO_INC}/CDF_DZero_Distribution.hh
IOHEADERS += ${IO_INC}/CDF_DZero_IOpoint.hh ${IO_INC}/CDF_DZero_IOfile.hh
IOHEADERS += ${IO_INC}/InfoForEfficiencyCalculation.hh

ALLSRCS := ${LIMITSRCS} ${IOSRCS}
ALLOBJS :=$(patsubst %.cpp,%.o,$(patsubst %.cxx,%.o,$(patsubst %.cc,%.o,${ALLSRCS})))

all: collie_io_bdt.exe collie_io_spectra.exe limit_calc.exe

collie_io_bdt.exe: collie_io_bdt.o
	gcc -o $@ $^ ${CXXFLAGS} ${LDFLAGS} ${LOADLIBS} ${EXAMPLELIBS} ${SEARCHDIRS} ${ROOTLIBS}

collie_io_spectra.exe: collie_io_spectra.o
	gcc -o $@ $^ ${CXXFLAGS} ${LDFLAGS} ${LOADLIBS} ${EXAMPLELIBS} ${SEARCHDIRS} ${ROOTLIBS}

limit_calc.exe: limit_calc.o
	gcc -o $@ $^ ${CXXFLAGS} ${LDFLAGS} ${LOADLIBS} ${EXAMPLELIBS} ${SEARCHDIRS} ${ROOTLIBS} -lCollieCLHEP

tidy::
	@${RM} core *.o *.exe

# target for removing all object files

clean:: tidy
	@${RM} make.depend make.depend.bak

# list of all source files

depend: make.depend

make.depend: Makefile
	touch make.depend
	makedepend -K -fmake.depend -- ${CXXFLAGS} -- ${ALLSRCS}  2> /dev/null

##include make.depend
# DO NOT DELETE
