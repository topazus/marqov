<table  align="center"><td align="center" width="9999">
<img src="doc/_static/logo_base.svg" align="center" width="120" alt="project icon">

<h1> MARQOV </h1>

<h4>A modern framework for classical spin models on general topologies</h4>

</td>
<tr>
<td align="left" width="9999" >
The MARQOV project is an extensible open-source framework for the massively parallel simulation and analysis of equilibrium lattice spin models. It enables users to study common problems in a unified framework but can also be molded to fit the needs for exploring advanced research problems.  The project website can be reached under [http://marqov.physik.uni-wuerzburg.de](http://marqov.physik.uni-wuerzburg.de/).<br><br>

For the physics motivation, Bramwel and Gingras put it best, when they said, that

>  Magnetic systems offer themselves as the
    ideal benchmark for generic concepts pertaining to collective phenomena in nature. This is
    due in part to the availability of a large variety
    of diverse magnetic materials that can be chosen to approximate simple theoretical “toy
    models” of collective behavior and, in part, to
    their ease of study by a battery of experimental
    techniques. <br>
    - Science 294(5546), 1495-1501 (2001).


MARQOV, together with its sibling pyMARQOV enables you to delve into this adventure and study the properties of these systems!


## Authors
* Manuel Schrauth (mschrauth@physik.uni-wuerzburg.de)
* Florian Goth (fgoth@physik.uni-wuerzburg.de)



## License
Every part of MARQOV is available under open licenses. Most parts are licensed under the GPLv3,
whereas some parts are even available under the more permissive MIT License.

## Platform/environment
MARQOV expects a common laptop/desktop like system and is routinely tested in a Linux environment.
Other POSIX-like platforms might also work.

### Tools/Dependencies
MARQOV is light on prerequisites. We require a C++14 compliant compiler and the cmake build system.
We use HDF5 for the output files, and will use the system HDF5 library if available.
If MPI is available, an MPI parallelized version of MARQOV will be built.
MARQOV features foreign language bindings that will be built if SWIG and the required language libraries are available.

## Configuration
MARQOV follows standard cmake procedures:

```console
cmake -E make_directory build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --target all --config Release
```

which will by default build all possible targets.
There is one switch available: -DDEBIASINTEGERS which can be used to build
a version of MARQOV that takes extra special care of eliminating bias in the generation of small random nubmers

## Testing
A small testsuite is available that can be accessed via 
- make test

## Usage
MARQOV is meant to be used as a library.
For a start have a look at the files in /demo.
There's mysimpleising and mysimpleheisenberg that show basic usage,
or have a look at our showcase MARQOVdemo.

Or just hack away to your heart's content by modifying main.cpp in src! There is a Makefile available, such that the compilation is particularly simply. The code can then be executed from any desired working directory and the output as well as log files are placed there (in respective subdirectories). We recommend using the base code folder (one level above the build directory) as the work directory. To compile and execute MARQOV's main program, use

```console
cd build 
make main
cd ..
./build/src/main
```

One may avoid switching through directories in the build process using make's -C flag

```console
make -C build/main
./build/src/main
```

</td>
</tr>
</table>

