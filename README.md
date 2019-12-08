# A library for integrating few-body motions
Please read this short manual before start to use the code. 
Any questions that are already answered here won't get reply from developers and GitHub Issue.

## User Guide
This c++11 library contains three components: Binary tree, AR and Hermite.
1. Binary tree: transformation between Kepler orbital parameters and positions and velocities of binary components; construction for hierarchical systems by using binary tree.
2. AR: algorithmical regularization or time-transformed explicit symplectic integrator with the slowdown method.
3. Hermite: a hybrid method combines the 4^{th} order Hermite integrator and the AR method.
The library can be used to create class objects embedded in other $N$-body code (see src directory) or used independently (see sample directory).
The complete user manual in html format is in docs/html (see GitHub page also) or pdf format in doct/pdf.

## Use as library
To use the code as a library in order to be implemented in other $N$-body code, developers can visit the directory src.
Three subdirectories exist:
1. Common: commonly used header files, including Float, List, ParticleGroup and BinaryTree.
1.1. Float.h: define the high-precision floating point type (Float)
1.2. List.h: define the template class List as a general container like std::vector
1.3. ParticleGroup.h: define a contain for a group of particles
1.4. BinaryTree.h: define a class to handle the transformation between Kepler orbital parameters and particle positions and velocities, also can be used to detect binaries and create hierarchical binary tree (also support hyperbolic orbits).
2. AR: algorithmical regularization or time-transformed explicit symplectic integrator with the slowdown method. The symplectic\_integrator.h is the header file to be included. In this file, two major classes are
2.1. SymplecticIntegrator: the main class that contain the particle data and integrator. 
2.2. SymplecticManager: the parameters used to control the integrator.
3. Hermite: the hybrid method that contain the Hermite integrator to integrating the global system and multiple AR integrators to integrate subsystems. The hermite.h is the header file to be included. In this file, two major classes are
3.1. HermiteIntegrator: the main class that contain the particle data and integrator. 
3.2. HermiteManager: the parameters used to control the integrator.

To use AR and Hermite, additional header files which define the interactions (calculation of the force and potential of pair interactions between singles and groups; calculation of perturbations for the slowdown method) should be provided separately, see sample codes for details.

## Use as independent code
In sample directories, three subdirectories provide sample codes using the library, which can be used as independent codes for integrating few-body motions.
1. Kepler: two tools for transformation between Kepler orbits and particle positions and velocities
1.1. keplerorbit: transformation between a group of binaries with mass, positions and velocities and a group of Kepler orbital parameters (semi-major axis, eccentricity, Euler angles and eccentric anomaly).
1.2. keplertree: generate a hierarchical Kepler binary tree from a group of particles (hyperbolic orbits can also be detected) or use a hierarchical orbital input to generate the group of particles (mass, positions and velocities).
2. AR: a $N$-body integrator using AR library. Four sample codes are generated: ar.logh, ar.logh.sd, ar.ttl, ar.ttl.sd. The suffix 'logh' and 'ttl' indicate the AR methods (see paper or docs for details); 'sd' indicates the slowdown method is included for inner binaries in a hierarchical systems.
3. Hermite: a $N$-body integrator combines Hermite and AR method, the code name is hermite.

In sample/input directory, users can find several input samples for using AR or Hermite codes.

## Questions and bug reports
For questions and bug reports, users can use the GitHub Issues to submit.
