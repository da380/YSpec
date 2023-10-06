MOD=../mod/
LOBJ=../obj/
BIN=../bin/
MV=mv -f
F95 = f95


FFLAGS = -J../mod  -O5




SRCS =  nrtype.f90 nrutil.f90 module_util.f90 module_spline.f90 \
	module_function.f90 module_fourier.f90 module_model.f90 \
	module_int.f90 \
	module_start.f90 module_solver.f90 


OBJS =  $(LOBJ)nrtype.o $(LOBJ)nrutil.o $(LOBJ)module_util.o \
	$(LOBJ)module_spline.o  $(LOBJ)module_function.o \
	$(LOBJ)module_fourier.o $(LOBJ)module_model.o  \
	$(LOBJ)module_int.o  $(LOBJ)module_start.o \
	 $(LOBJ)module_solver.o



MODS =  $(MOD)nrtype.mod $(MOD)nrutil.mod $(MOD)module_util.mod \
	$(MOD)module_spline.mod $(MOD)module_function.mod \
	$(MOD)module_fourier.mod $(MOD)module_model.mod \
	$(MOD)module_int.mod  $(MOD)module_start.f90 \
	$(MOD)module_solver.mod   




all:    yspec      \
        yspec_pre  \
        yspec_pro  \
        yspec_cal  \
        yspec_test \
        int_test1  \
        int_test2

clean: 
	rm ../bin/yspec
	rm ../bin/yspec_pre
	rm ../bin/yspec_cal
	rm ../bin/yspec_pro
	rm ../bin/yspec_test 
	rm ../bin/int_test1
	rm ../bin/int_test2
	rm ../obj/*
	rm ../mod/*

yspec: $(BIN)yspec $(OBJS)
$(BIN)yspec: $(OBJS) $(LOBJ)yspec.o
	$(F95) $(FFLAGS) $(LOBJ)yspec.o $(OBJS) -o $(BIN)yspec

$(LOBJ)yspec.o: yspec.f90 
	$(F95) $(FFLAGS) yspec.f90 -c
	$(MV) yspec.o $(LOBJ)


yspec_pre: $(BIN)yspec_pre $(OBJS)
$(BIN)yspec_pre: $(OBJS) $(LOBJ)yspec_pre.o
	$(F95) $(FFLAGS) $(LOBJ)yspec_pre.o $(OBJS) -o $(BIN)yspec_pre

$(LOBJ)yspec_pre.o: yspec_pre.f90 
	$(F95) $(FFLAGS) yspec_pre.f90 -c
	$(MV) yspec_pre.o $(LOBJ)


yspec_pro: $(BIN)yspec_pro $(OBJS)
$(BIN)yspec_pro: $(OBJS) $(LOBJ)yspec_pro.o
	$(F95) $(FFLAGS) $(LOBJ)yspec_pro.o $(OBJS) -o $(BIN)yspec_pro

$(LOBJ)yspec_pro.o: yspec_pro.f90 
	$(F95) $(FFLAGS) yspec_pro.f90 -c
	$(MV) yspec_pro.o $(LOBJ)




yspec_cal: $(BIN)yspec_cal $(OBJS)
$(BIN)yspec_cal: $(OBJS) $(LOBJ)yspec_cal.o
	$(F95) $(FFLAGS) $(LOBJ)yspec_cal.o $(OBJS) -o $(BIN)yspec_cal

$(LOBJ)yspec_cal.o: yspec_cal.f90 
	$(F95) $(FFLAGS) yspec_cal.f90 -c
	$(MV) yspec_cal.o $(LOBJ)



yspec_test: $(BIN)yspec_test $(OBJS)
$(BIN)yspec_test: $(OBJS) $(LOBJ)yspec_test.o
	$(F95) $(FFLAGS) $(LOBJ)yspec_test.o $(OBJS) -o $(BIN)yspec_test

$(LOBJ)yspec_test.o: yspec_test.f90 
	$(F95) $(FFLAGS) yspec_test.f90 -c
	$(MV) yspec_test.o $(LOBJ)


int_test1: $(BIN)int_test1 $(OBJS)
$(BIN)int_test1: $(OBJS) $(LOBJ)int_test1.o
	$(F95) $(FFLAGS) $(LOBJ)int_test1.o $(OBJS) -o $(BIN)int_test1

$(LOBJ)int_test1.o: int_test1.f90 
	$(F95) $(FFLAGS) int_test1.f90 -c
	$(MV) int_test1.o $(LOBJ)


int_test2: $(BIN)int_test2 $(OBJS)
$(BIN)int_test2: $(OBJS) $(LOBJ)int_test2.o
	$(F95) $(FFLAGS) $(LOBJ)int_test2.o $(OBJS) -o $(BIN)int_test2

$(LOBJ)int_test2.o: int_test2.f90 
	$(F95) $(FFLAGS) int_test2.f90 -c
	$(MV) int_test2.o $(LOBJ)



$(OBJS): $(SRCS)	
	$(F95)   $(FFLAGS) -c $(SRCS)
	$(MV)  *.o $(LOBJ)


