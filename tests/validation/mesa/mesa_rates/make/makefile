# Builds tests/validation/mesa/mesa_rates/mesa_rates against $(MESA_DIR) (same MESA as the bbq references).
ifndef MESA_DIR
$(error MESA_DIR enviroment variable is not set)
endif
include $(MESA_DIR)/utils/makefile_header

TEST_DIR = ..
TEST_SRC_DIR = ../src
LOAD_OTHER = -L$(MESA_LIB_DIR) -lnet -leos -lrates -lchem -lutils $(LOAD_MESA_NUMERICS)
OBJS = mesa_rates.o

$(TEST_DIR)/mesa_rates : $(OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(OBJS) $(LOAD_OTHER) $(LD_FLAGS)

%.o: $(TEST_SRC_DIR)/%.f90
	$(TEST_COMPILE) $(FCfree) $<

clean:
	-@rm -f *.o *.mod *.smod ../mesa_rates
