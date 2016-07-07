PYTHON  = python3
SWIGDIR = swig
IDDIR   = idcpp

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
	-rm -rf .depend *.o *.so dist build idcpp __pycache__ idcpp.egg-info $(SWIGDIR)/*_wrap.cxx $(SWIGDIR)/idcpp.py
