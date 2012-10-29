include ../../Make.inc
LIBDIR=../../lib
PSBLIBDIR=$(PSBLASDIR)/lib
PSBINCDIR=$(PSBLASDIR)/include
MLDLIBDIR=$(MLDDIR)/lib
MLDINCDIR=$(MLDDIR)/include
HERE=..

FINCLUDES=$(FMFLAG).. $(FMFLAG)$(LIBDIR) $(FMFLAG)$(MLDINCDIR)  $(FMFLAG)$(PSBINCDIR) 


F90OBJS=psb_d_cp_dsc_from_coo.o \
psb_d_cp_dsc_from_fmt.o \
psb_d_cp_dsc_to_coo.o \
psb_d_cp_dsc_to_fmt.o \
psb_d_dsc_aclsum.o \
psb_d_dsc_allocate_mnnz.o \
psb_d_dsc_arwsum.o \
psb_d_dsc_colsum.o \
psb_d_dsc_cp_from.o \
psb_d_dsc_csgetblk.o \
psb_d_dsc_csgetptn.o \
psb_d_dsc_csgetrow.o \
psb_d_dsc_csmm.o \
psb_d_dsc_csmv.o \
psb_d_dsc_csnm1.o \
psb_d_dsc_csnmi.o \
psb_d_dsc_csput.o \
psb_d_dsc_cssm.o \
psb_d_dsc_cssv.o \
psb_d_dsc_get_diag.o \
psb_d_dsc_mold.o \
psb_d_dsc_mv_from.o \
psb_d_dsc_print.o \
psb_d_dsc_reallocate_nz.o \
psb_d_dsc_reinit.o \
psb_d_dsc_rowsum.o \
psb_d_dsc_scal.o \
psb_d_dsc_scals.o \
psb_d_dsc_trim.o \
psb_d_mv_dsc_from_coo.o \
psb_d_mv_dsc_from_fmt.o \
psb_d_mv_dsc_to_coo.o \
psb_d_mv_dsc_to_fmt.o \
mld_dainvk_impl.o \
mld_d_ainvk_solver_apply.o \
mld_d_ainvk_solver_apply_vect.o \
mld_dainvt_impl.o \
mld_d_ainvt_solver_apply.o \
mld_d_ainvt_solver_apply_vect.o \
mld_d_aorth_impl.o \
mld_d_aorth_solver_apply.o \
mld_d_aorth_solver_apply_vect.o \
mld_dinvk_copyin.o \
mld_dinvk_copyout.o \
mld_dinvk.o \
mld_dinvt_copyin.o \
mld_dinvt_copyout.o \
mld_dinvt.o \
mld_dsparse_ainvk.o \
mld_dsparse_ainvt.o

#mld_dainvt_impl.o mld_dainvk_impl.o mld_d_aorth_impl.o \
# psb_d_dsc_impl.o


OBJS=$(F90OBJS) 

LIBNAME=libainv_prec.a

lib: $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)


veryclean: clean

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

