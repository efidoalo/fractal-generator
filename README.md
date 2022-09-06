# fractal-generator
A GUI program (tested on ubuntu20.04) that generates and displays different fractals depending on variables set by the user.
The GNU MPFR library (https://mpfr.loria.fr/#intro)  is used to dynamically select significand bit length for calcuating the mandelbrot set with high precision and zoom in capabilities. The program also uses my own vector data structure implementation.
Compile and link with the following command: gcc `pkg-config --cflags gtk+-3.0` -o fractals fractals.c `pkg-config --libs gtk+-3.0` -lm ~/Documents/Containers/C/vector.o -lmpfr -lgmp
