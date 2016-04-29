FC=ifort

FFLAGS = -c -save -w -V -WB -opt-prefetch -double-size 64 -assume buffered_io -ip -O3 -xhost -mkl 
LDFLAGS = -zero

synthe:
	$(FC) $(FFLAGS)  -o atlas7lib.o  atlas7v.for
	$(FC) $(FFLAGS) -o xnfpelsyn.o xnfpelsyn.for
	$(FC) $(LDFLAGS) -o xnfpelsyn.exe xnfpelsyn.o atlas7lib.o
	$(FC) $(FFLAGS)  -o synbeg.o synbeg.for
	$(FC) $(LDFLAGS)  -o synbeg.exe synbeg.o
	$(FC) $(FFLAGS) -o rmolecasc.o rmolecasc.for
	$(FC) $(LDFLAGS)  -o rmolecasc.exe rmolecasc.o
	$(FC) $(FFLAGS)  -o rgfall.o  rgfall.for
	$(FC) $(LDFLAGS)  -o rgfall.exe  rgfall.o
	$(FC) $(FFLAGS) -o rpunchbin.o rpunchbin.for
	$(FC) $(LDFLAGS) -o rpunchbin.exe rpunchbin.o
	$(FC) $(FFLAGS)  -o synthe.o synthe.for
	$(FC) $(LDFLAGS) -o synthe.exe synthe.o
	$(FC) $(FFLAGS) -o spectrv.o spectrv.for
	$(FC) $(LDFLAGS) -o spectrv.exe spectrv.o atlas7lib.o
	$(FC) $(FFLAGS) -o broadenx.o broadenx.for
	$(FC) $(LDFLAGS) -o broadenx.exe broadenx.o
	$(FC) $(FFLAGS) -o rotate.o rotate.for
	$(FC) $(LDFLAGS) -o rotate.exe rotate.o
	$(FC) $(FFLAGS) -o syntoascanga.o syntoascanga.for
	$(FC) $(LDFLAGS)  -o syntoascanga.exe syntoascanga.o

testing:
	$(FC) $(FFLAGS)  -o synthe.o synthe.for
	$(FC) $(LDFLAGS) -o synthe.exe synthe.o

install:
	/bin/mv *.exe ../bin

clean:
	rm -f *.o
