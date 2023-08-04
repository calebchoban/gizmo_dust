FLAGS = -DWSS_CIE_COOL -DCHEMCOOL -DCHEMISTRYNETWORK=5 -g #-r8 -traceback

FC = gfortran $(FLAGS)

OBJS  = coolinmo.o cheminmo.o spline.o cool_func.o photoinit_ism.o \
        dvode.o evolve_abundances_dust.o rate_eq_simple.o jac.o cool_util.o \
        const_rates.o validate_rates.o calc_photo.o calc_temp.o \
	compute_gamma.o wss_z_collis.o harness.o 
#calc_shield_test.o 
#test_harness.o

test: $(OBJS)
	$(FC) -o test $(OBJS)

clean:
	rm -f $(OBJS)
