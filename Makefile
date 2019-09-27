# python bin
PYBIN =/home/gzn/anaconda3/bin/python

default:
	export PYTHONPATH=`pwd`/lib/python
	$(PYBIN) setup.py install 
	@echo "#!"${PYBIN} -u > cak.py; cat _CALYPSO_ANALYSIS.py >> cak.py
	@chmod a+x cak.py
	@if [ -d lib64 ]; then echo export PYTHONPATH=`pwd`/lib64/python':$$PYTHONPATH' > caly.sh; \
	   	else echo export PYTHONPATH=`pwd`/lib/python':$$PYTHONPATH' > caly.sh;fi
	@echo export PATH=`pwd`':$$PATH' >> caly.sh
	@echo
	@echo please add \" source `pwd`/caly.sh \" in $$HOME/.bashrc
	@echo
	

.PHONY : clean

clean:
	-rm -rf build lib lib64 caly.sh

