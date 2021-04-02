# This is an commentary line in a makefile
# Start of the makefile


#
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996- .
#
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

SIESTA_ARCH=x86_64-ifort12-mpiifort

FPP=
FPP_OUTPUT=
FC=mpiifort
FC_SERIAL=ifort
RANLIB=ranlib

SYS=nag

SP_KIND=4
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)

FFLAGS=-I/share/apps/intel/mkl/include -I/share/apps/intel/mkl/include/intel64/lp64 -g -O2
FPPFLAGS= -DFC_HAVE_FLUSH -DFC_HAVE_ABORT -DMPI -DGRID_DP -DPHI_GRID_SP
LDFLAGS=

ARFLAGS_EXTRA=

FCFLAGS_fixed_f=
FCFLAGS_free_f90=
FPPFLAGS_fixed_F=
FPPFLAGS_free_F90=

BLAS_LIBS=-L/share/apps/intel/mkl/lib/intel64 -lmkl_solver_lp64 -lmkl_intel_lp64 -lpthread -limf -lm
LAPACK_LIBS=-L/share/apps/intel/mkl/lib/intel64 -lmkl_sequential -lmkl_core
BLACS_LIBS=-L/share/apps/intel/mkl/lib/intel64 -lmkl_blacs_intelmpi_lp64
SCALAPACK_LIBS=-L/share/apps/intel/mkl/lib/intel64 -lmkl_scalapack_lp64

COMP_LIBS=

NETCDF_LIBS=
NETCDF_INTERFACE=

LIBS=$(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(NETCDF_LIBS)

#SIESTA needs an F90 interface to MPI
#This will give you SIESTA's own implementation
#If your compiler vendor offers an alternative, you may change
#to it here.
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/share/apps/intel/impi/4.0.3.008/intel64/include

 

IF90 = ifort
OBJS =  $(OUT_DIR)/ConstantsUnit.o \
        $(OUT_DIR)/RootsUnit.o \
        $(OUT_DIR)/1RandUtils.o \
        $(OUT_DIR)/0Integration.o \
        $(OUT_DIR)/0Integration2.o \
        $(OUT_DIR)/0IntialParameters.o \
        $(OUT_DIR)/2SubroutineFunctions.o \
        $(OUT_DIR)/3BasicVariablesAndMethods.o \
        $(OUT_DIR)/4EmittersBlackBody.o \
        $(OUT_DIR)/Chandra_Diffuse_Reflection.o \
        $(OUT_DIR)/5ScatterPhotons.o \
        $(OUT_DIR)/6ScatterDistance.o \
        $(OUT_DIR)/7EstimationUnit.o \
        $(OUT_DIR)/8photons.o \
        $(OUT_DIR)/9SubFuncsForDiffuseReflection.o \
        $(OUT_DIR)/main.o
 
        #$(OUT_DIR)/5CrossSection.o \
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#$(OUT_DIR)/RootsUnit.o: $(OUT_DIR)/ConstantsUnit.o
#$(OUT_DIR)/ellfunctionsUnit.o: $(OUT_DIR)/RootsUnit.o
#$(OUT_DIR)/BlcoordinateUnit.o: $(OUT_DIR)/ellfunctionsUnit.o
$(OUT_DIR)/0IntialParameters.o: $(OUT_DIR)/ConstantsUnit.o
$(OUT_DIR)/0Integration.o: $(OUT_DIR)/ConstantsUnit.o
$(OUT_DIR)/0Integration2.o: $(OUT_DIR)/ConstantsUnit.o
$(OUT_DIR)/2SubroutineFunctions.o: $(OUT_DIR)/ConstantsUnit.o $(OUT_DIR)/1RandUtils.o
$(OUT_DIR)/3BasicVariablesAndMethods.o: $(OUT_DIR)/0Integration.o $(OUT_DIR)/2SubroutineFunctions.o \
             $(OUT_DIR)/0IntialParameters.o
$(OUT_DIR)/4EmittersBlackBody.o: $(OUT_DIR)/3BasicVariablesAndMethods.o
$(OUT_DIR)/Chandra_Diffuse_Reflection.o: $(OUT_DIR)/3BasicVariablesAndMethods.o
$(OUT_DIR)/6ScatterPhotons.o: $(OUT_DIR)/Chandra_Diffuse_Reflection.o
$(OUT_DIR)/6ScatterDistance.o: $(OUT_DIR)/6ScatterPhotons.o
$(OUT_DIR)/7Photons.o: $(OUT_DIR)/6ScatterDistance.o $(OUT_DIR)/4EmittersBlackBody.o 
$(OUT_DIR)/8SubFunc_ThomScat_Diffuse_Reflection.o: $(OUT_DIR)/6ScatterDistance.o \
			 $(OUT_DIR)/7Photons.o $(OUT_DIR)/Chandra_Diffuse_Reflection.o  
$(OUT_DIR)/main.o: $(OUT_DIR)/8SubFunc_ThomScat_Diffuse_Reflection.o 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

$(shell if [ ! -e $(OUT_DIR) ];then mkdir -p $(OUT_DIR); fi)
$(shell if [ ! -e $(MOD_DIR) ];then mkdir -p $(MOD_DIR); fi) 


$(bin): $(OBJS)
	$(FC) $^ -o ../$@ $(LDFLAGS) $(LDLIBS) 
 

$(MOD_DIR)/%.mod: %.o %.f90


#$(OUT_DIR)/%.o: $(YNOGK_DIR)/%.f90 Makefile
#	$(FC) -c -CB -O2 $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90) -I$(OUT_DIR) -module $(MOD_DIR) $< -o $@

$(OUT_DIR)/%.o: $(LEMON_DIR)/%.f90 Makefile
	$(FC) -c -CB -O2 $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90) -I$(OUT_DIR) -module $(MOD_DIR) $< -o $@
 
$(OUT_DIR)/%.o: %.f90
	$(FC) -c -CB -O2 $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90) -I$(OUT_DIR) -module $(MOD_DIR) $< -o $@

#$(bin): $(OBJS) $(LIB_DIR)/libynogk.a
#$(bin): $(OBJS) $(LIB_DIR)/libynogk.a
#	$(FC) $^ -o ../$@ $(LDFLAGS) $(LDLIBS) -I$(YNOGK_DIR)

#$(YNOGK_DIR)/%.o: %.f90
#	$(FC) -c -CB -O2 $<

#%.mod: %.o %.f90
#%.o: %.f90 $(LIB_DIR)/libynogk.a
#%.o: %.f90 $(LIB_DIR)/libynogk.a
#	$(FC) -c -CB -O2 $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90) $< -I$(YNOGK_DIR)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# static library
#$(LIB_DIR)/libynogk.a:
#	$(MAKE) -C $(YNOGK_DIR) $@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
cleanall:
	rm $(bin) $(OBJ) *.mod


#${BIN}: ${OBJS}
#	${IF90} -o  ${BIN} ${OBJS} 

#%.mod: %.o %.f90  
#%.o: %.f90   
#	${IF90} -c -O3 -check bounds $<
 
#clean:
#	rm *.o *.mod compton
# End of the makefile
# End of the makefile
