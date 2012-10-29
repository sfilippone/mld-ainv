include ../../Make.inc
LIBDIR=../../lib
PSBLIBDIR=$(PSBLASDIR)/lib
PSBINCDIR=$(PSBLASDIR)/include
MLDLIBDIR=$(MLDDIR)/lib
MLDINCDIR=$(MLDDIR)/include
HERE=..

FINCLUDES=$(FMFLAG).. $(FMFLAG)$(LIBDIR) $(FMFLAG)$(MLDINCDIR)  $(FMFLAG)$(PSBINCDIR) 


F90OBJS=mld_dainvt_impl.o mld_dainvk_impl.o mld_d_aorth_impl.o  psb_d_dsc_impl.o


OBJS=$(F90OBJS) 

LIBNAME=libainv_prec.a

lib: $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)


veryclean: clean

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

