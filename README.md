Secuential implementation of Finite Diferences for heat equation solution in 1D and 2D.
---------------------------------------------------------------------------------------
for now  we just code explicit solution for heat equation in 1D and 2.


###Table of contents
**[Finite Diferences for heat equation in 1D](#1D)**

**[Finite Diferences for heat equation in 2D](#2D)**

##Finite Diferences for heat equation in 1D.

Explicit solution in 1D directory

### Dependencies
* [Armadillo](http://arma.sourceforge.net/ )
* [octave](https://www.gnu.org/software/octave/ )

### Compile and run

    mkdir build
    cd build
    cmake ../
    make
    ./fdtd_1D < ../in

keep in mind that **fdtd_1D** executable receives the following parameters that you can put in file **in**:

* xa: init position in X
* xb: last position in X
* ta: init value in time
* tb: last value in time
* N: number of point to calculate heat value
* T: number of steps in time.
* gamma: heat difusion value in equation.

### Runing example
remove build directory and exec in terminal

    ./run.sh

the output example image is "image.jpg"

### example Docs
this program had been tested with these papers:

[heat Equation example 1](http://sistemas.fciencias.unam.mx/~gcontreras/joomla15/images/stories/teoria.pdf)

[heat Equation example 2](http://www.icmc.usp.br/CMS/Arquivos/arquivos_enviados/BIBLIOTECA_113_RT_327.pdf)

##Finite Diferences for heat equation in 2D (rectangular objects)

Explicit solution in 2D directory

### Dependecies
* [Armadillo](http://arma.sourceforge.net/ )
* [Gnuplot](http://www.gnuplot.info/)

###Compile and Run
in 2D directory type:

    ./run.sh


executable file is **./build/fdtd_2d** and its input file is **in2d**
and have de following parameters:

* xa: init position in X
* xb: last postition in X
* ya: init position in y
* yb: last position in y
* ta: init time
* tb: last time
* Nx: discretization points in X
* Ny: discretization points in Y
* T: discretization points in time
* k: heat difusion value in equation

###Output
output values are init state and last state of heat in object.

###Runing and ploting example

in 2D directory delete build dir (if you have it ) and run:

    ./run.sh
    ./plot_Example

this is a little example of rectangular object with a center heat source and plot in files **time0.png** and **timeT.png** the init and end heat state of object


###References
[two dimensional heat equation with FD](http://geodynamics.usc.edu/~becker/teaching/557/problem_sets/problem_set_fd_2dheat.pdf)
