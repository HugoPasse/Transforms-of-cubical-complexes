all:
	cythonize -a -i embedded_cubical_complex.pyx
	mv embedded_cubical_complex.cpython-38-x86_64-linux-gnu.so embedded_cubical_complex.so

clean:
	rm -rf build
	rm embedded_cubical_complex.so
	rm embedded_cubical_complex.cpp