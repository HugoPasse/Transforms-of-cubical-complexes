# Transforms-of-cubical-complexes

Code in this repository implements algorithms to compute hybrid transforms, Radon transform and Euler characteristic transform on cubical complexes. You have two choices, you can use the file `Embedded_cubical_complex.h` as a header file in a C++ program, or you can compile the python module to use it in a python program.

## Python module compilation

To use our C++ code in python we use cython. To compile the module you will need the  [`Gudhi`](https://gudhi.inria.fr/) library headers. Once you've downloaded them, modify the file called `setup.py` :
```python
libs = ["/your_path_to_gudhi"]
```
 Then you can compile the module by typing :
```
python3 setup.py build_ext --inplace
```

Don't forget to add `~/your_module.so` to your python path (use `export PYTHONPATH=~/your_module.so`) or to put the `your_module.so` file into you python's module directory.
