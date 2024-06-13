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


A summary of the input parameters is as follows:

-  prefix for output files – a string used for labelling the outputs of the program. For exam-
ple, if this string is yspec.out, then the output files will be yspec.out.1, yspec.out.2, . . . ,
yspec.out.n, where n denotes the number of receivers, and the ordering of these files matches the
input order of the receivers (see below).


2. earth model – name of the file listing the earth model to be used in the calculations. The input
model must be in the ‘deck format’ used in the program minos. The model can have a fluid inner
core, and the program takes this into account without any explicit input.
3. attenuation switch – set to one to include attenuation in the calculations, and to zero if atten-
uation should not be included. Note that the form of attenuation is the logarithmic dispersion
relation used in the parameterization of the PREM model. Furthermore, the code implements this
attenuation exactly, and does not depend upon first-order perturbation theory (as is, usually, done
in the case of normal mode summation).
4. gravitation – The code can incorporate the effects of gravitation in a number of ways. Set param-
eter to zero to ignore gravitation completely, equal to one to employ the Cowling approximation,
and to two to enable full self-gravitation. Using either no gravitation or the Cowling approximation
will be a little faster, but not substantially so.
5. output – Determines the form of the output seismograms. Set to zero to get displacements, equal
to one for velocities, and equal to two for accelerations. In all cases the outputs are in SI units.
6. potential and tilt corrections – set equal to one to make potential and tilt corrections to the
spherioidal solutions, and zero otherwise. Note, if the previous option must be set equal to two (i.e.
to output accelerations) for these corrections to be applied.
7. lmin – minimum value of spherical harmonic degree l used in the calculations – probably should
always be set equal to one.
8. lmax – Maximum spherical harmonic degree for the calculations. The appropriate choice of this
parameter depends on the maximum frequency and also source-receiver geometry. If lmax is too
small the seismograms will be ringy prior to the arival of the direct p-wave – if in doubt make this
parameter bigger and see if anything changes.
9. fmin – minimum frequency for the calculations in mHz. Must be greater than zero – the equations
are singular in the zero-frequency limit, and this code cannot calculate static displacements.
10. fmax – maximum frequency for teh calculations in mHz.
