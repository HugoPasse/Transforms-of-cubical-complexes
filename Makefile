all:
	cythonize -a -i embedded_cubical_complex.pyx
	mv embedded_cubical_complex.*.so embedded_cubical_complex.so

clean:
	rm -rf build
	rm embedded_cubical_complex.so
	rm embedded_cubical_complex.cpp
	rm embedded_cubical_complex.html
