# Chris's Experiment Scripts - Synthetic FT data analysis - Makefile
# Written by Christopher Thomas.

MATLAB=matlab
MATFLAGS=-nodisplay -nodesktop -nosplash
MATGUIFLAGS=-nodesktop -nosplash

default:
	@echo "Targets:   run rungui force gallery clean"


force: clean run gallery


clean:
	rm -f plots/*

checkdirs:
	./checkdirs.sh


run: checkdirs
	nice -n +10 $(MATLAB) $(MATFLAGS) \
		-r "dbstop if error; run('do_test_ft.m'); exit;"

rungui: checkdirs
	nice -n +10 $(MATLAB) $(MATGUIFLAGS) \
		-r "dbstop if error; run('do_test_ft.m'); exit;"

gallery: gallery4

gallery6:
	sync
	cd plots; makegallery.pl --width=16% --nofolders *png; cd ..

gallery5:
	sync
	cd plots; makegallery.pl --width=19% --nofolders *png; cd ..

gallery4:
	sync
	cd plots; makegallery.pl --width=24% --nofolders *png; cd ..

gallery3:
	sync
	cd plots; makegallery.pl --width=32% --nofolders *png; cd ..

gallery2:
	sync
	cd plots; makegallery.pl --width=49% --nofolders *png; cd ..

#
# This is the end of the file.
