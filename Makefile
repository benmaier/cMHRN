clean:
	-rm -f *.o
	-rm -f $(TARGET)

clean_all:
	make clean
	sudo make pyclean
	make matclean

pyclean:
	-rm -f *.so
	-rm -rf *.egg-info*
	-rm -rf ./tmp/
	-rm -rf ./build/

matclean:
	-rm -rf ./matlabbuild/

py:
	sudo python setup.py develop
