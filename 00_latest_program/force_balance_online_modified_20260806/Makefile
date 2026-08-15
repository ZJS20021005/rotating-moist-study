#Makefile for parallel compiling  
FC=h5pfc -r8 -check bounds -traceback#-mkl#-fpp -check bounds
FC+=-O3 -fma -fpp -qno-openmp #-ip -ipo
FC += -I$(HOME)/lib/fftw-3.3.7/include

# flow config
#FC+=-DMOVIE
FC+=-DPLANEMOVIE
FC+=-DSTRUCMOVIE

# Library
LINKS=-lhdf5_fortran -lhdf5
LINKS+=-L$(HOME)/lib/fftw-3.3.7/lib -lfftw3 #-lfftw3_omp
LINKS+=-L/home/linux/eva/lib/fftw-3.3.7/lib -lfftw3
#LINKS+=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LINKS+=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#LINKS+=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

#============================================================================ 
# make PROGRAM   
#============================================================================
OBJECTS_RAW=inimov_hdf_ycut.o mkmov_hdf_ycut.o inimov_hdf_xcut.o mkmov_hdf_xcut.o inimov_hdf_zcut.o mkmov_hdf_zcut.o mkmov_mean.o mkmov_field.o mkmov_dsal.o mkmov_qtens.o mkmov_qvap.o \
			initParams.o cfl.o cflr.o coetar.o cordin.o continua.o divg.o divgck.o \
			divgloc.o gcurv.o compute_qsat_from_b.o hdnlq1.o hdnlq2.o hdnlq3.o hdnlsa.o hdnlqvap.o \
            indic.o inimov.o inirea.o initia.o inqpr.o invtrq1.o invtrq2.o \
            invtrq3.o invtrsa.o invtrqvap.o matrix_transpose.o meshes.o \
            mpi_routines.o openfi.o param.o papero.o phcalc.o phini.o \
            prcalc.o solq1k.o solq2k.o solq3k.o solsak.o solqvap.o solxi_endwall.o solxj_endwall.o solxi_fs.o solxj_fs.o solxi_ns.o solxj_ns.o solxri.o solxrj.o\
            solxi_periodic.o solxj_periodic.o solxri_periodic.o solxrj_periodic.o \
            tsch.o updvp.o avgvar.o force_balance_online.o vmaxv.o nusse_walls.o caldiffusivity.o StatRoutines.o \
            aggregation_scale_diagnostics.o \
            mse_vortex_diagnostics.o \
			coetarr.o indicr.o write_field.o velbc.o\
			mgrd_idc.o mgrd_mem.o mgrd_velitp.o mgrd_dsalc.o mgrd_qvapc.o \

PROGRAM=simexec
MODULES=obj/param.o
MODLOCS=mod/
OBJLOCS=obj/
OBJECTS=$(addprefix $(OBJLOCS),$(OBJECTS_RAW))
#============================================================================ 
# Linking    
#============================================================================
$(PROGRAM):objlocs modlocs $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) $(LINKS) -module $(MODLOCS) -o $@

#============================================================================
#  Dependencies
#============================================================================
obj/param.o: ./param.f90
	$(FC) -c -o $@ -module $(MODLOCS) $<

obj/tsch.o: ./tsch.f90
	$(FC) -c -o $@ -module $(MODLOCS) $<

obj/gcurv.o: ./gcurv.f90
	$(FC) -c -o $@ -module $(MODLOCS) $<

obj/papero.o: ./papero.f90
	$(FC) -c -o $@ -module $(MODLOCS) $<

obj/%.o: particle_ibm/%.f90 $(MODULES)
	$(FC) -c -o $@ -module $(MODLOCS) $<

obj/%.o: fluid_solver/%.f90 $(MODULES)
	$(FC) -c -o $@ -module $(MODLOCS) $<

obj/%.o: output/%.f90 $(MODULES)
	$(FC) -c -o $@ -module $(MODLOCS) $<
#============================================================================
#  Clean up
#============================================================================
clean :
	rm obj/*.o 
	rm mod/*.mod
	rm $(PROGRAM)

#============================================================================
#  Definitions
#============================================================================
.PHONY: objlocs
objlocs: ${OBJLOCS}
${OBJLOCS}:
	mkdir -p ${OBJLOCS}
.PHONY: modlocs
modlocs: ${MODLOCS}
${MODLOCS}:
	mkdir -p ${MODLOCS}
