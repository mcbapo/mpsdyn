# mpsdyn
Working version of MPS library (variational GS and excited state search, time evolution and a few models)

There are two folders: 
    src (containing all the code) 
    doc, with some (automatically generated via doxygen) documentation. To start browsing, open doc/html/index.html, you can browse the structure of the code, classes, etc.

The first thing you need to do is to install the primme static library libprimme.a (https://github.com/primme/primme) in src/libs. 
    Depending on your architecture, it should be in src/libs/linux or src/libs/macos 

If primme is not available, the code can run using arpack, but this has had some issues in the past and is less efficient.
    
It also requires a BLAS and LAPACK installation, which you probably have.
For early tests and some old programs, also the arpack library is used (this is to be removed at some point). 
The path to these libraries should be added to the variable LINKDIRS in Makefile.arch under the suitable option (corresponding to the architecture). If possible, using MKL libraries is recommended. 

Right now, the code runs, as it is provided, in a Linux cluster (Scientific Linux) and some Mac computers, but one may need to touch the Makefile.arch (in src/bin) for a particular architecture. 

Once this is done, in src/bin, one should copy the file Makefile.newuser to Makefile.user.[your_user_name] (in src/bin/  changing the last part of the name to your actual username). Then you can compile and run some tests (Please, start by compiling "`make test1`" and running "`./test1`" to test basic tensor functions. If it runs without failing, the tensor library, and interface to primme etc are working)

The interesting programs are in src/programs/  and different subdirectories

Compilation instructions should be in the Makefile.user.[your_user_name] file. For instance

make ising

compiles the program programs/ising/testIsing.cpp
This computes the GS of the Ising model, variationally  and with imaginary time evolution. It requires a list of parameters:

    ./ising L J g h D M delta outputFile append

- `L`(system size) 
- `J g h` (Hamiltonian parameters) 
- `D` (bond dimension) 
- `M` (number of time steps for imag. t evol.) 
- `delta` (time step) 
- `outputFile` (name, including relative path, of output file) 
- `append` (0/1 whether the file is to be recreated or appended to)


To run other programs, you will need to provide a configuration file with a set of properties (required ones are specified in the .cpp file in programs ). As an example, there is some in folder config/
