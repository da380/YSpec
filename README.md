# YSpec

The program yspec calculates synthetic seismograms in spherically symmetric earth models using the
direct radial integration method introduced by Friedrich & Dalkolmo (1995), and extended in Al-Attar
& Woodhouse (2008) to incorporate the effects of self-gravitation.

## Installation

The code has no dependencies, and requires only a Fortran90 compiler. 

Clone the repository:

`git clone https://github.com/da380/YSpec.git`

Configure the build:

`cmake -S YSpec -B YSpecBuild -DCMAKE_BUILD_TYPE=Release`

Compile the code:

`cmake --build YSpecBuild`

## Running the code


The program is run using the command:

`yspec input.in`

where `input.in` is the name of the parameter file. And example of such a file can be found within the `examples` folder within the build directory. 


