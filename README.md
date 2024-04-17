# SDAR: A library for solving few-body problems
Please read this short manual before start to use the code. 
Any questions that are already answered here won't get reply from developers and GitHub Issue.
You can find the [SDAR documentation](https://lwang-astro.github.io/SDAR/docs/html/index.html) written in the Doxygen style.

# About the version
SDAR code is maintained across multiple branches, and for users seeking stability, it is advisable to obtain the released version (see tags on GitHub webpage). The master branch undergoes regular updates to introduce new features and address bugs. Users can opt for this version if they find the new features beneficial.

Other branches may not function correctly, and users are advised against interacting with them unless they have consulted the developer and comprehended the implementations and existing issues.

Users can retrieve the code version in VERSION file.
The major modes of versions are as follows:

- **Release mode (r)**: This mode signifies a stable version recommended for general use. Release versions are infrequently updated, with the initial release currently being prepared.
- **Master mode**: This mode is suitable for most conditions and is expected to function correctly. Users can use it for production runs.
- **Develop mode (e)**: In this mode, the code should operate correctly but may contain bugs that lead to unexpected results. It is advisable to consult developers before use.
- **Test mode (test)**: This mode is unverified for proper functioning and should not be employed for production purposes.

## User Guide
This c++11 library contains three components: _BinaryTree_, _AR_ and _Hermite_.
1. _BinaryTree_: transformation between Kepler orbital parameters and positions and velocities of binary components; construction for hierarchical systems by using binary tree.
2. _AR_ algorithmical regularization or time-transformed explicit symplectic integrator with the slowdown method.
3. _Hermite_: a hybrid method combines the 4^{th} order Hermite integrator and the AR method.
The library can be used to create class objects embedded in other _N_-body code (see src directory) or used independently (see sample directory).
The complete user manual in html format is in docs/html (see GitHub page also) or pdf format in doct/pdf.

## Use as library
To use the code as a library in order to be implemented in other _N_-body code, developers can visit the directory src.
Three subdirectories exist:
1. Common: commonly used header files, including Float, List, ParticleGroup and BinaryTree.
    1. Float.h: define the high-precision floating point type (Float)
    2. List.h: define the template class List as a general container like std::vector
    3. ParticleGroup.h: define a contain for a group of particles
    4. BinaryTree.h: define a class to handle the transformation between Kepler orbital parameters and particle positions and velocities, also can be used to detect binaries and create hierarchical binary tree (also support hyperbolic orbits).
2. AR: algorithmical regularization or time-transformed explicit symplectic integrator with the slowdown method. The symplectic\_integrator.h is the header file to be included. In this file, two major classes are
    1. SymplecticIntegrator: the main class that contain the particle data and integrator. 
    2. SymplecticManager: the parameters used to control the integrator.
3. Hermite: the hybrid method that contain the Hermite integrator to integrating the global system and multiple AR integrators to integrate subsystems. The hermite.h is the header file to be included. In this file, two major classes are
    1. HermiteIntegrator: the main class that contain the particle data and integrator. 
    2. HermiteManager: the parameters used to control the integrator.

To use AR and Hermite, additional header files which define the interactions (calculation of the force and potential of pair interactions between singles and groups; calculation of perturbations for the slowdown method) should be provided separately, see sample codes for details.

## Use as independent code
In sample directories, three subdirectories provide sample codes using the library, which can be used as independent codes for integrating few-body motions.
Use `make` to compile the sample codes and `make install` to install the excutable fieles.
The default install path is the bin directory in user's home path. 
To change the install path, please modify `INSTALL_PATH` in Makefile in each sample subdirectory.

Here are the description of three subdirectories:
1. Kepler: two tools for transformation between Kepler orbits and particle positions and velocities.
    1. keplerorbit: transformation between a group of binaries with mass, positions and velocities and a group of Kepler orbital parameters (semi-major axis, eccentricity, Euler angles and eccentric anomaly).
    2. keplertree: generate a hierarchical Kepler binary tree from a group of particles (hyperbolic orbits can also be detected) or use a hierarchical orbital input to generate the group of particles (mass, positions and velocities).
2. AR: a _N_-body integrator using AR library. Six sample executable files are generated: ar.logh, ar.logh.sd.a ar.log.sd.t, ar.ttl, ar.ttl.sd.a ar.ttl.sd.t. The meaning of suffixes:
    1. 'logh' and 'ttl' indicate the AR methods (see paper or docs for details).
    2. 'sd' indicates the slowdown method is included for inner binaries in a hierarchical systems. 
    3. 'a' means that the slowdown method is only used for inner binaries in a hierarchical system; 't' means that the slowdown method is used hierachically for inner and outer binaries in a system.
3. Hermite: a _N_-body integrator combines Hermite and AR method, the code name is hermite.

In sample/input directory, users can find several input samples for using AR or Hermite codes.

## Questions and bug reports
For questions and bug reports, users can use the GitHub Issues to submit.
