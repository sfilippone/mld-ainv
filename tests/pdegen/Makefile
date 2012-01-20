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

PDOBJS=ppde.o  data_input.o 
PSOBJS=spde.o  data_input.o 
MLOBJS=mldpde.o data_input.o
EXEDIR=./runs

all: mldpde
# ppde spde
mldpde: $(MLOBJS)
	$(F90LINK) $(MLOBJS) -o mldpde $(AINV_LIB) $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv mldpde $(EXEDIR)

ppde: $(PDOBJS)
	$(F90LINK) $(PDOBJS) -o ppde $(AINV_LIB) $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS) -Wl,-pthread
	/bin/mv ppde $(EXEDIR)


spde: $(PSOBJS)
	$(F90LINK) $(PSOBJS) -o spde $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv spde $(EXEDIR)

ppde.o spde.o: data_input.o 


clean: 
	/bin/rm -f $(PDOBJS) $(EXEDIR)/ppde $(PSOBJS) $(EXEDIR)/spde 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



