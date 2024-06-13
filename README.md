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

-  `prefix for output files` – a string used for labelling the outputs of the program. For exam-
ple, if this string is `yspec.out`, then the output files will be `yspec.out.1`, `yspec.out.2`, . . . ,
`yspec.out.n`, where `n` denotes the number of receivers, and the ordering of these files matches the
input order of the receivers (see below).


- `earth model` – name of the file listing the earth model to be used in the calculations. The input
model must be in the "deck format" used in the program minos. The model can have a fluid inner
core, and the program takes this into account without any explicit input.

- `attenuation switch` – set to `1` to include attenuation in the calculations, and to `0` if atten-
uation should not be included. Note that the form of attenuation is the logarithmic dispersion
relation used in the parameterization of the PREM model. Furthermore, the code implements this
attenuation exactly, and does not depend upon first-order perturbation theory (as is, usually, done
in the case of normal mode summation).

- `gravitation` – The code can incorporate the effects of gravitation in a number of ways. Set param-
eter to `0` to ignore gravitation completely, equal to `1` to employ the Cowling approximation,
and to `2`  to enable full self-gravitation. Using either no gravitation or the Cowling approximation
will be a little faster, but not substantially so. The method in the non-gravitating and self-gravitating
cases have been carefully benchmarked, but this has not been done for the Cowling approximation.

- `output` – Determines the form of the output seismograms. Set to `0` to get displacements, equal
to `1` for velocities, and equal to `2` for accelerations. In all cases the outputs are in SI units.

- `potential and tilt corrections` – set equal to `1` to make potential and tilt corrections to the
spherioidal solutions, and `0` otherwise. Note, if the previous option must be set equal to `2` (i.e.
to output accelerations) for these corrections to be applied.

- `lmin` – minimum value of spherical harmonic degree `l` used in the calculations – probably should
always be set equal to `0`.

- `lmax` – Maximum spherical harmonic degree for the calculations. The appropriate choice of this
parameter depends on the maximum frequency and also source-receiver geometry. If `lmax` is too
small the seismograms will be ringy prior to the arival of the direct p-wave – if in doubt make this
parameter bigger and see if anything changes. As a rule of thumb for PREM or similar models, `lmax ~ 15 fmax`
with `fmax` the maximum frequency in mHz.

- `fmin` – minimum frequency for the calculations in mHz. Must be greater than zero – the equations
are singular in the zero-frequency limit, and this code cannot calculate static displacements.

- `fmax` – maximum frequency for the calculations in mHz.

- `length of time series` – length of the desired time series in minutes as measured from the origin
time of the event.

- `time step` – the constant time-step in seconds desired for the output seismograms.
  
- `f11, f12, f21, f22` – parameters for cosine high and low pass filters applied prior to calculation
of the inverse Fourier transform. A high-pass filter is applied between `f11` and `f12` , and a low-pass
filter between `f21` and `f22` . If all these parameters are set to `0`, no filtering will be done.
Note, however, that some form of low-pass filtering must be done to remover the effects of the sharp
cut-off of the calculations at the frequency `fmax`.

- `source depth, latitude and longitude` – location of point source (depth in km, angles in de-
grees).

- `Mrr , . . . , Mpp` – moment tensor components $M_{rr} , . . . , M_{φφ}$ in Nm refered to the standard spherical
polar co-ordinate system (r, θ, φ) defined relative to the pole of the model. So, locally the `r` direction
is vertically upwards, the `θ` direction is south, and the `φ` direction is east.

16. receiver depth – depth of all the different receivers in km. Note that the receiver can be placed
below the source.
17. number of receivers – the total number of different receiver positions. Can be as large as you
like, but if too big the program can run out of memory and will crash. Computations for multiple
receivers involves only slightly more work than those for one receiver, so having lots of different
receivers is fine.
18. receiver latitudes and longitudes – list of the different receiver latitudes and longitudes in
the form
lat1 lon1
lat2 lon2
lat3 lon3
with all angles being in degrees. The ordering of the output files corresponds to that of this list.
