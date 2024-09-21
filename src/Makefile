default: CIMI

CIMI:
	make EXE

CIMI_NO_FLC:
	make EXE_NO_FLC


OBJECT=\
       	plasmasphere.f90 \
       	w2k.f T96_01.f TS04c.f \
       	geopack_2008.f \
       	Lstar2.f90 \
	rcm6.f90 \
	gmresm.f90 

ifort=ifort -r8 -132 -save -O3 -mcmodel=medium 
#ifort=ifort -r8 -132 -save -O3 -qopenmp
debug=ifort -r8 -132 -save -O3 -traceback -g -check -debug -mcmodel=medium
gfort=gfortran -fdefault-real-8 -ffixed-line-length-132 -fno-automatic -fno-align-commons -O4 


FLC_TEST:
	@cp ModCurvScatt.f90 ./src_unit_test/
	@cd src_unit_test;	make EXE_FLC
	@cd src_unit_test;	rm ModCurvScatt.f90
	cp src_unit_test/unit_test_curv_radius.exe ./run_test/
	cp src_unit_test/unit_test_curv_scatt.exe ./run_test/
	cp tools/plot_d_curv_scatt.pro ./run_test/
	@cd run_test;	./unit_test_curv_radius.exe
	@cd run_test;	./unit_test_curv_scatt.exe

EXE:
	${ifort} -o cimi.exe ModCurvScatt.f90 cimi.f90 ${OBJECT}
	#${ifort} -o cimi.exe ModCurvScatt.f90 cimi_main.f90 ${OBJECT}

debug:
	${debug} -o cimi.exe ModCurvScatt.f90 cimi.f90 ${OBJECT}

EXE_NO_FLC:
#	${gfort} -o cimi_no_flc.exe ${OBJECT} 
	${ifort} -o cimi_no_flc.exe cimi_no_flc.f90 ${OBJECT} 

clean:
	rm -rf *.mod
	rm *.exe
