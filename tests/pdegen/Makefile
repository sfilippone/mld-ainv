AINVDIR=../..
include $(AINVDIR)/Make.inc
PSBLIBDIR=$(PSBLASDIR)/lib/
PSBINCDIR=$(PSBLASDIR)/include
MLDINCDIR=$(MLDDIR)/include
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIB=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
AINVLIBDIR=$(AINVDIR)/lib
AINV_LIB=-L$(AINVLIBDIR) -lainv_prec
PSBLAS_LIB= -L$(PSBLIBDIR)  -lpsb_util  -lpsb_base
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AINVLIBDIR) $(FMFLAG)$(MLDINCDIR)  $(FMFLAG)$(PSBINCDIR) $(FIFLAG).

D3DOBJS=ppde3d.o data_input.o  pde3d_exp_mod.o pde3d_gauss_mod.o pde3d_base_mod.o
EXEDIR=./runs

all: ppde3d
# ppde spde
ppde3d: $(D3DOBJS)
	$(F90LINK) $(D3DOBJS) -o ppde3d $(AINV_LIB) $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde3d $(EXEDIR)

ppde3d.o:  data_input.o  pde3d_exp_mod.o pde3d_gauss_mod.o pde3d_base_mod.o

clean: 
	/bin/rm -f $(PDOBJS) $(EXEDIR)/ppde $(PSOBJS) $(EXEDIR)/spde 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



