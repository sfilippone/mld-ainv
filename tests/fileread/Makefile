MLDAINVDIR=../..
include $(MLDAINVDIR)/Make.inc
PSBLIBDIR=$(PSBLASDIR)/lib/
PSBINCDIR=$(PSBLASDIR)/include
MLDMODDIR=$(MLDDIR)/modules
MLDINCDIR=$(MLDDIR)/include
MLDLIBDIR=$(MLDDIR)/lib
PSBMODDIR=$(PSBLASDIR)/modules
MLD_LIB=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
MLDAINVINCDIR=$(MLDAINVDIR)/include
MLDAINVMODDIR=$(MLDAINVDIR)/modules
MLDAINVLIBDIR=$(MLDAINVDIR)/lib
MLDAINV_LIB=-L$(MLDAINVLIBDIR) -lainv_prec
PSBLAS_LIB= -L$(PSBLIBDIR)  -lpsb_util  -lpsb_base
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDAINVMODDIR)  $(FMFLAG)$(MLDAINVINCDIR)  $(FMFLAG)$(MLDINCDIR) $(FMFLAG)$(PSBINCDIR) $(FIFLAG).  $(FMFLAG)$(MLDMODDIR)  $(FMFLAG)$(PSBMODDIR) 

DFSOBJS=df_sample.o data_input.o
SFSOBJS=sf_sample.o data_input.o
CFSOBJS=cf_sample.o data_input.o
ZFSOBJS=zf_sample.o data_input.o

EXEDIR=./runs

all: df_sample
#sf_sample cf_sample zf_sample

df_sample: $(DFSOBJS)
	$(FLINK) $(LINKOPT) $(DFSOBJS) -o df_sample \
	$(MLDAINV_LIB) $(MLD_LIB)  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv df_sample $(EXEDIR) 

sf_sample: $(SFSOBJS)
	$(FLINK) $(LINKOPT) $(SFSOBJS) -o sf_sample \
	$(MLDAINV_LIB) $(MLD_LIB)  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv sf_sample $(EXEDIR) 

cf_sample: $(CFSOBJS)
	$(FLINK) $(LINKOPT) $(CFSOBJS) -o cf_sample \
	$(MLDAINV_LIB) $(MLD_LIB)  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv cf_sample $(EXEDIR) 

zf_sample: $(ZFSOBJS)
	$(FLINK) $(LINKOPT) $(ZFSOBJS) -o zf_sample \
	$(MLDAINV_LIB) $(MLD_LIB)  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv zf_sample $(EXEDIR) 

sf_sample.o: data_input.o
df_sample.o: data_input.o
cf_sample.o: data_input.o
zf_sample.o: data_input.o

clean: 
	/bin/rm -f $(DFOBJS) $(ZFOBJS) $(SFSOBJS)  $(DFSOBJS) \
        *$(.mod)  $(EXEDIR)/df_bench $(EXEDIR)/zf_bench \
	$(EXEDIR)/sf_sample $(EXEDIR)/cf_sample  \
	$(EXEDIR)/df_sample $(EXEDIR)/zf_sample

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

