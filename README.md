# YSpec

The program yspec calculates synthetic seismograms in spherically symmetric earth models using the
direct radial integration method introduced by Friedrich & Dalkolmo (1995), and extended in Al-Attar
& Woodhouse (2008) to incorporate the effects of self-gravitation.

## Installation

Clone the repository:

`git clone https://github.com/da380/YSpec.git`

Configure the build:

`cmake -S YSpec -B YSpecBuild -DCMAKE_BUILD_TYPE=Release`

Compile the code:

`cmake --build YSpecBuild`


