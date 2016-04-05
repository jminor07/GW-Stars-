# a makefile
.SUFFIXES:
.SUFFIXES: .out .o .c .e .r .F .f .y .l .s .p .vmsf .for .vmsc .com
#       Teach make how to handle Fortran programs with C-preprocessor lines
.F.o: ; f77 $(FFLAGS) -c $*.F
#       Teach make how to make files for export to a VMS system
.vmsf.for: ; cp $*.vmsf $*.for
.vmsc.com: ; cp $*.vmsc $*.com
.F.for: ; /lib/cpp $(VMSFLAGS) $*.F > $*.for
.f.for: ; cp $*.f $*.for

padask: padask.f
	f77 padask.f -o padask /space/www/dial/code/libchem.a \
        /usr/a/theory/mods/one/code_fitsio/libfitsio.sol.a  -lm
wwwask: wwwask.f
	f77 wwwask.f -o wwwask /space/www/dial/code/libchem.a \
        /usr/a/theory/mods/one/code_fitsio/libfitsio.sol.a  -lm
wwwpad: wwwpad.f
	f77 wwwpad.f -o wwwpad exportpad2/code/libinterp.sol.a \
        /usr/a/theory/mods/one/code_fitsio/libfitsio.sol.a  -lm
driver: driver.f
	f77 driver.f -o driver
