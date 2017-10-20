MLDAINVDIR=../..
include $(MLDAINVDIR)/Make.inc
PSBLIBDIR=$(PSBLASDIR)/lib/
PSBINCDIR=$(PSBLASDIR)/include
PSBMODDIR=$(PSBLASDIR)/modules
MLDINCDIR=$(MLDDIR)/include
MLDLIBDIR=$(MLDDIR)/lib
MLDMODDIR=$(MLDDIR)/modules
MLD_LIB=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
MLDAINVINCDIR=$(MLDAINVDIR)/include
MLDAINVMODDIR=$(MLDAINVDIR)/modules
MLDAINVLIBDIR=$(MLDAINVDIR)/lib
MLDAINV_LIB=-L$(MLDAINVLIBDIR) -lainv_prec 
PSBLAS_LIB= -L$(PSBLIBDIR)  -lpsb_util  -lpsb_base
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDAINVMODDIR)  $(FMFLAG)$(MLDAINVINCDIR) $(FMFLAG)$(MLDINCDIR)  $(FMFLAG)$(PSBINCDIR) $(FIFLAG). $(FMFLAG)$(MLDMODDIR)  $(FMFLAG)$(PSBMODDIR) 


D3DOBJS=ppde3d.o data_input.o  pde3d_exp_mod.o pde3d_gauss_mod.o pde3d_base_mod.o \
	genpde_mod.o
D2DOBJS=ppde2d.o data_input.o  pde2d_exp_mod.o pde2d_gauss_mod.o pde2d_base_mod.o \
	genpde_mod.o

EXEDIR=./runs

all: ppde3d ppde2d

ppde3d: $(D3DOBJS)
	$(FLINK) $(D3DOBJS) -o ppde3d $(MLDAINV_LIB) $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde3d $(EXEDIR)

ppde3d.o:  data_input.o  pde3d_exp_mod.o pde3d_gauss_mod.o pde3d_base_mod.o genpde_mod.o

ppde2d: $(D2DOBJS)
	$(FLINK) $(D2DOBJS) -o ppde2d $(MLDAINV_LIB) $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde2d $(EXEDIR)

ppde2d.o:  data_input.o  pde2d_exp_mod.o pde2d_gauss_mod.o pde2d_base_mod.o genpde_mod.o

clean: 
	/bin/rm -f $(D2DOBJS) $(D3DOBJS) *.$(.mod)

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



