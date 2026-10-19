# YSpec

The program yspec calculates synthetic seismograms in spherically symmetric earth models using the
direct radial integration method introduced by Friedrich & Dalkolmo (1995), and extended in Al-Attar
& Woodhouse (2008) to incorporate the effects of self-gravitation.

## Installation

The code requires a Fortran compiler (e.g. gfortran) and the [FFTW3](https://www.fftw.org) library
(e.g. the `libfftw3-dev` package), which is used for the transform to the time domain. If CMake does
not find FFTW3, set `FFTW3_INCLUDE_DIR` and `FFTW3_LIBRARY`. If the compiler supports OpenMP, the
loop over frequencies is run in parallel; set `-DYSPEC_USE_OPENMP=OFF` to build a serial version. The
results are identical in either case.

Clone the repository:

`git clone https://github.com/da380/YSpec.git`

Configure the build:

`cmake -S YSpec -B YSpecBuild -DCMAKE_BUILD_TYPE=Release`

Compile the code:

`cmake --build YSpecBuild`

Run the regression tests (optional):

`ctest --test-dir YSpecBuild`

The tests compare the outputs for a few small cases in the `tests` folder with stored references.

## Running the code


The program is run using the command:

`yspec input.in`

where `input.in` is the name of the parameter file. An example of such a file can be found within the `examples` folder within the build directory. 

The number of threads can be set using the environment variable `OMP_NUM_THREADS`.

The parameter file is read by counting lines, so its layout must match the example exactly: two
header lines, and then each value preceded by exactly one blank line and one comment line. File
names that contain a `/` must be enclosed in quotes (e.g. `"../models/prem.200"`), as otherwise
they are truncated at the first `/`.


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
spheroidal and radial solutions, and `0` otherwise. Note that the previous option must be set equal
to `2` (i.e. to output accelerations) for these corrections to be applied.

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
grees). Latitudes of the source and receivers are taken to be geographic, and are converted to
geocentric latitudes using a flattening of about 1/297. A point lying exactly at either pole is
moved 0.02 degrees away from it so that its azimuths are defined.

- `Mrr , . . . , Mpp` – moment tensor components $M_{rr} , . . . , M_{\phi}$ in Nm refered to the standard spherical
polar co-ordinate system (r, θ, φ) defined relative to the pole of the model. So, locally the $r$ direction
is vertically upwards, the $\theta$ direction is south, and the $\phi$ direction is east. The source time
function is a Heaviside step, i.e. the moment tensor switches on at the origin time and stays on.

- `receiver depth` – depth of all the different receivers in km. Note that the receiver can be placed
below the source. This depth cannot lie within a fluid part of the model.

- `number of receivers` – the total number of different receiver positions. Can be as large as you
like, but if too big the program can run out of memory and will crash. Computations for multiple
receivers involves only slightly more work than those for one receiver, so having lots of different
receivers is fine.

- `receiver latitudes and longitudes` – list of the different receiver latitudes and longitudes in
the form

`lat1 lon1`

with all angles being in degrees. The ordering of the output files corresponds to that of this list.

## Outputs

For each receiver is output an ascii file containing the three desired three component seismogram. The
format is:

`time   vertical-component   north-component   east-component`

where all outputs are in SI units, and the form of the seismogram (displacement, velocity, or acceler-
ation) is determined by parameter 5 described above.

## Changes

The following errors have been corrected. Seismograms computed with earlier versions of the code
will differ, mostly in the horizontal components:

- The toroidal contribution to the horizontal component along the source-receiver great circle used
  $X_{lm}$ in place of $X_{lm}/\sin\theta$. This component was wrong by an amount that vanishes at an
  epicentral distance of 90 degrees but grows towards the source and the antipode (e.g. by about 25%
  of its amplitude at 20 degrees, and by about the size of the signal at 5 degrees).
- Toroidal modes with $l = 1$ were left out. They change the seismograms by 0.1–1%.
- The spheroidal source terms left out the imaginary (attenuating) part of the elastic moduli at the
  source. The error is of order $1/Q$ at the source depth: about 0.03% for a crustal source in PREM, and
  about 1% for a source in the low-Q zone.
- The free-air correction was not applied to the radial ($l = 0$) modes when potential and tilt
  corrections were requested.
- A receiver at exactly 90 degrees north was placed at the south pole.
- Several minor issues: degree-to-radian constants were truncated, a few quantities were computed in
  single precision, some arrays were used before being set, and `lmax = 0` wrote out of bounds.
- The constant $\pi$ was only accurate to single precision. Through the Fourier transform this caused
  errors of about $10^{-6}$ of the peak amplitude.
- The spherical Bessel functions used for the starting values of the integrations have been replaced
  by a single, accurate continued fraction, and the depth at which the integrations start is now
  chosen without reference to PREM-specific frequency thresholds and more conservatively. The
  seismograms now differ from those obtained by starting the integrations at the base of the model
  by about $10^{-8}$ of their peak values.

The code has also been tidied up. In particular, the routines that were taken from Numerical Recipes,
whose source code cannot be redistributed, have been replaced, with FFTW now used for the Fourier
transform, and unused code has been removed. Error conditions now end the program with a non-zero
exit status, and the output time series run from 0 to the requested length. The remaining
differences from the outputs of earlier versions are no larger than the accuracy of the numerical
integrations.




