default: CIMI

rundir_path='./run/'

CIMI:
	if test -d $(rundir_path); then\
		echo "the run directory exist";\
	else \
		make rundir; \
	fi
	@cd src; make
	cp src/cimi.exe ./run/

CIMI_NO_FLC:
	make EXE_NO_FLC
	cp cimi_no_flc.exe ./run/


rundir:
	rm -rf run
	mkdir run
	cp inputs/* ./run
	rm ./run/D_*.dat 
	cd run; ln -s ../inputs/D_* ./
	cd run; ln -s ./cimi.dat.4ion ./cimi.dat
	cp tools/*.pro ./run
	#inputs/* tools/*.pro run

debug:
	@cd src; make debug

clean:
	rm -rf run
	rm src/*.mod
	rm src/cimi.exe
