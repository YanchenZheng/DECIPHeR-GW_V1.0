PROGRAM        = DECIPHeR_v1_UK_res.exe
PROGRAM_i        = DECIPHeR_v1_i_UK_res.exe

.PHONY: $(PROGRAM) $(PROGRAM_i) clean cleanall

all: $(PROGRAM)

$(PROGRAM):
	cd debug && make
	cp debug/$(PROGRAM) .

$(PROGRAM_i):
	cd intel && make
	cp intel/$(PROGRAM) ./$(PROGRAM_i)

clean: /;
	-rm $(PROGRAM) $(PROGRAM_i)
	-rm *.mod *.o
	cd debug && make clean
	cd intel && make clean
