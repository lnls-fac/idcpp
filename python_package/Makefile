PYTHON  = python3
SWIGDIR = swig
IDDIR   = insertion_devices

.NOTPARALLEL:

install:
	rm -f $(IDDIR)/VERSION
	cp $(SWIGDIR)/__init__.py $(IDDIR)
	ln -s -r VERSION $(IDDIR)
	$(PYTHON) setup.py install

develop:
	rm -f $(IDDIR)/VERSION
	cp $(SWIGDIR)/__init__.py $(IDDIR)
	ln -s -r VERSION $(IDDIR)
	$(PYTHON) setup.py develop

clean:
	-rm -rf .depend *.o *.so dist build insertion_devices __pycache__ insertion_devices.egg-info $(SWIGDIR)/*_wrap.cxx $(SWIGDIR)/insertion_devices.py
