include ../Make.inc
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules
PSBLIBDIR=$(PSBLASDIR)/lib
PSBINCDIR=$(PSBLASDIR)/include
PSBMODDIR=$(PSBLASDIR)/modules
MLDLIBDIR=$(MLDDIR)/lib
MLDINCDIR=$(MLDDIR)/include
MLDMODDIR=$(MLDDIR)/modules
HERE=.

FINCLUDES=$(FMFLAG). $(FMFLAG)$(LIBDIR) $(FMFLAG)$(MLDINCDIR)  $(FMFLAG)$(PSBINCDIR) \
 $(FMFLAG)$(MLDMODDIR)  $(FMFLAG)$(PSBMODDIR) 


MODOBJS=mld_d_invt_solver.o mld_base_ainv_mod.o \
	mld_d_invk_solver.o mld_d_ainv_solver.o \
	mld_d_biconjg_mod.o mld_ainv_mod.o mld_d_base_ainv_mod.o mld_d_ainv_tools_mod.o

#F90OBJS=mld_dainvt_impl.o mld_dainvk_impl.o  psb_d_dsc_impl.o mld_d_aorth_bld_mod.o 


OBJS=$(F90OBJS) $(MODOBJS) 

LIBMOD=
LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=libainv_prec.a


libcopy: lib impld
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) *$(.mod) $(MODDIR)

lib: $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

impld: lib
	$(MAKE) -C impl LIBNAME=$(LIBNAME)

mld_d_base_ainv_mod.o mld_d_biconjg_mod.o: mld_base_ainv_mod.o mld_d_ainv_tools_mod.o
mld_d_ainv_solver.o mld_d_invt_solver.o mld_d_invk_solver.o: mld_d_base_ainv_mod.o
mld_ainv_mod.o: mld_d_ainv_solver.o mld_d_invt_solver.o mld_d_invk_solver.o 
veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: iclean
	/bin/rm -f $(OBJS) $(LOCAL_MODS) *$(.mod)
iclean:
	$(MAKE) -C impl clean
