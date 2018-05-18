default: 
	make python

clean:
	-rm -f *.o
	make pyclean

clean_all:
	make clean
	make pyclean

pyclean:
	-rm -f *.so
	-rm -rf *.egg-info*
	-rm -rf ./tmp/
	-rm -rf ./build/

python:
	pip install -e ../cMHRN --no-binary :all:

grootinstall:
	source /home/bfmaier/.bash_profile
	/opt/python36/bin/pip3.6 install --user ../cMHRN

groot:
	git fetch
	git pull
	make clean
	make grootinstall
