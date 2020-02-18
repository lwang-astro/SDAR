/*! \mainpage Introduction
  
\section Algorithm

The slow-down time-transformed explicit symplectic method (SDAR) combines the benefit of the symplectic integrator which conserves the Hamiltonian and angular momentum and the high efficiency of the slow-down method to handle the long-term evolution of hierarchical systems and close encounters.
The details of the algorithm can be found in Wang, Nitadori & Makino (2020).
Here the idea is briefly introduced.

The slow-down method is introduced by <A HREF="http://adsabs.harvard.edu/abs/1996CeMDA..64..197M">Mikkola & Aarseth (1996)</A>.
For a perturbed binary, by artificially slowing down the orbital motion (scaling the time) while keeping the orbital parameters unchanged, the external perturbation is effectively enlarged.
In such case, the effect of perturbation on one binary orbit can represent the average effect of several orbits.
When the perturbation is weak, this method can properly approximate the secular evolution.

The symplectic integrator can conserve the Hamiltonian and the angular momentum of a system.
Thus it is very suitable for simulating the long-term evolution of a system.
However, it requires a constant integration step.
In the classical symplectic method, time step, \f$ \mathrm{d} t \f$, is also the integration step, \f$ \mathrm{d} s \f$.
This leads to a low efficiency in integrating an eccentric Kepler orbit.
In order to be accurately enough, \f$ \mathrm{d} t \f$ has to be fixed to the smallest value determined at the pericenter.
The solution is to apply a time transformation, \f$\mathrm{d} t=g \mathrm{d} s\f$, which decouples \f$\mathrm{d} s\f$ and \f$\mathrm{d} t\f$ with a function \f$g\f$ (see <A HREF="http://www.sciencedirect.com/science/article/pii/S0168927497000615"> Hairer, 1997 </A>).
Thus, \f$ \mathrm{d} t \f$ can vary to avoid the issue of low efficiency while \f$ \mathrm{d} s \f$ is fixed to keep the symplectic property.
This method can be described by the extended phase space Hamiltonian:

(1) \f$ \Gamma(\mathbf{P},\mathbf{Q}) = H'(\mathbf{P},\mathbf{Q}) \frac{d t}{d s} = g(\mathbf{Q},\mathbf{P}) (H(\mathbf{p},\mathbf{q},t) + Pt) \f$

where \f$ H(\mathbf{p},\mathbf{q},t) \f$ is the standard Hamiltonian.
The extended phase-space vector, \f$ (\mathbf{P},\mathbf{Q}) \f$, is defined by the standard coordinate and momentum pair \f$(\mathbf{p},\mathbf{q})\f$ with an additional pair of the coordinate, \f$ t \f$, and the conjugate momentum, \f$ Pt = - H(\mathbf{p},\mathbf{q},0) \f$ (negative initial Hamiltonian).

The equation of motion with the differential on \f$ s \f$ describe the symplectic map in the extended phase-space.
The disadvantage is that the time transformation usually results in an inseparable Hamiltonian and only the expensive implicit integrator can be used.
<A HREF="http://adsabs.harvard.edu/abs/1999MNRAS.310..745M">Mikkola & Tanikawa (1999)</A> and <A HREF="http://adsabs.harvard.edu/abs/1999AJ....118.2532P">Preto & Tremaine (1999)</A> find a solution by defining a specific type of time transformation function:

(2) \f$ g(\mathbf{Q},\mathbf{P}) = \frac{f(T(\mathbf{P})) - f(-U(\mathbf{Q}))}{T(\mathbf{P}) + U(\mathbf{Q})} \f$

that the Hamiltonian can be written in a separable style, in order to use the explicit symplectic integrator.
Especially, for an isolated binary system, if \f$ f(x) = log(x) \f$, so that $\mathrm{d} t \approx \mathrm{d} s/ |U|$, where $|U|$ is the absolute value of the binary potential energy (similar like the Burdet-Heggie time transformation), the integrator behaviors dramatically well for the Kepler orbit, i.e., the numerical trajectory follows the exact one with a phase error of time.

\subsection H_sec Hamiltonian in Extended Phase Space

Defining the general coordinates as \f$ \mathbf{q} = \{q_i\}, (i=1,n) \f$ with freedom of \f$n\f$ and corresponding general momentums ad \f$ \mathbf{p} \f$,  The Hamiltonian equations is:

(1) \f$ \frac{d \mathbf{q}}{d t} = \frac{\partial H}{\partial \mathbf{p}}\f$; \f$ \frac{d \mathbf{p}}{d t} = - \frac{\partial H}{\partial \mathbf{q}} \f$

Here the dt is used as a differetial varaible. 
For the propuse as we discussed above, we want to use a new variable \f$s\f$ replacing the function of time \f$t\f$. 
In this case, the time is treated as a new general coordinate. 
And the corresponding time momentum \f$Pt\f$ should be also added.

We extend the coordiantes to \f$ \mathbf{Q} = (t, \mathbf{q}) \f$ and the momentums to \f$ \mathbf{P} = (Pt, \mathbf{p})\f$  with total freedom of \f$2(n+1)\f$.

The new Hamiltonian \f$H'\f$ should also satisfy the Hamiltonian equations (1). 
Especially for \f$(t, Pt)\f$, we can get:

(2) \f$ \frac{d t}{d t} = \frac{\partial H'}{\partial Pt} = 1 \f$; \f$ \frac{d Pt}{d t} = - \frac{\partial H'}{\partial t} = - \frac{\partial H}{\partial t}\f$

From first equation of (2), we find the \f$H'\f$ linearly depend on \f$Pt\f$, thus \f$H'\f$ can be the form as \f$ H' = H + Pt \f$.
The second equation indicates that the time evolution of \f$Pt\f$ is equal to the negative energy change of the system.
Thus the value of \f$Pt\f$ at the time \f$t\f$ can be \f$-H(t)\f$.

We want to write Hamiltonian equations with new differetial variable \f$ ds\f$.
Defining \f$ g(\mathbf{Q},\mathbf{P}) = \frac{dt}{ds} \f$, we can rewrite (1) with \f$ds\f$ and extended coordinates \f$(\mathbf{Q}, \mathbf{P})\f$ as:

(3) \f$ \frac{d \mathbf{Q}}{d s} = g(\mathbf{Q},\mathbf{P}) \frac{\partial H'}{\partial \mathbf{P}} \f$;   \f$ \frac{d \mathbf{P}}{d s} = - g(\mathbf{Q},\mathbf{P}) \frac{\partial H'}{\partial \mathbf{Q}} \f$

However, we need to have the Hamiltonian equations the same form as original, thus we need to find another Hamiltonian \f$\Gamma(\mathbf{P},\mathbf{Q})\f$ that satisfy the Hamiltonian equations:

(4) \f$ \frac{d \mathbf{Q}}{d s} = \frac{\partial \Gamma}{\partial \mathbf{P}} \f$; \f$ \frac{d \mathbf{P}}{d s} = -\frac{\partial \Gamma}{\partial \mathbf{Q}} \f$

To find correct \f$\Gamma(\mathbf{P},\mathbf{Q})\f$, we go back to the Principle of least action which is used to derive the Lagrangian equations.
The relation between (standard) Hamiltonian \f$H(\mathbf{p},\mathbf{q},t)\f$ and Lagrangian \f$L(\mathbf{p},\mathbf{q},t)\f$ is 

(5) \f$ H(\mathbf{p},\mathbf{q},t) = \sum_{i=1}^n p_i \dot{q_i} - L(\mathbf{p},\mathbf{q},t) \f$

The Principle of least action require the action

(6) \f$ S = \int_{t_1}^{t_2} L(\mathbf{p},\mathbf{q},t) dt = \int_{t_1}^{t_2} \left[ \sum_{i=1}^n p_i \dot{q_i} - H(\mathbf{p},\mathbf{q},t) \right] dt \f$ 

should take the mimimum path, thus any function variation \f$ \delta S \f$ should makes \f$ S + \delta S\f$ increase.
Thus when \f$ \delta L(\mathbf{p},\mathbf{q},t) = 0 \f$, this condition is satisfied. This leads to the Lagrangian equations and also the Hamitonian equations.

Here the integration takes from \f$ t_1 \f$ to \f$ t_2 \f$ and the time is used as integration variable. 
Now we treat \f$(t, Pt)\f$ as new coordinate and momemtum, \f$H'\f$ as new Hamitonian, and use \f$s\f$ as new integration variable.
Then \f$S\f$ can be rewrited as:

(7) \f$ S = \int_{s_1}^{s_2} \left[ \sum_{i=1}^n p_i \frac{d q_i} {d s} + Pt \frac{d t}{d s} - (H(\mathbf{p},\mathbf{q},t) + Pt) \frac{d t}{d s} \right] ds = \int_{s_1}^{s_2} \left[ \sum_{i=1}^{n+1} P_i \frac{d Q_i}{d s} - H'(\mathbf{P},\mathbf{Q}) \frac{d t}{d s}\right] ds \f$

It is obvious that when

(8) \f$ \Gamma(\mathbf{P},\mathbf{Q}) = H'(\mathbf{P},\mathbf{Q}) \frac{d t}{d s} = g(\mathbf{Q},\mathbf{P}) (H(\mathbf{p},\mathbf{q},t) + Pt) \f$

The formula (7) become the same form as (6). Then with Principle of least action, the Hamiltonian equation (4) can be derived.
We call the \f$ \Gamma(\mathbf{P},\mathbf{Q}) \f$ is the Hamiltonian in the extended phase space \f$ (\mathbf{P},\mathbf{Q}) \f$

The Hamiltonian in extended phase space \f$\Gamma\f$ is also useful for analyzing the systems where Hamiltonian \f$H\f$ explicitly depends on time and is not conserved. 
Since time become a coordinate in \f$\Gamma\f$, \f$\frac{\partial \Gamma}{\partial s}\f$ is zero thus \f$ \Gamma\f$ become conserved quantity. 
The method dealing with closed system can be used with Hamiltonian in extended phase space.

\subsection T_sec Time transformation for Separable Hamiltonian

With the Hamiltonian in extended phase space, we can integrate the equation of motions with step \f$ ds \f$ by choosing a kind of \f$g(\mathbf{Q},\mathbf{P})\f$. If we choose a \f$g(\mathbf{Q},\mathbf{P})\f$ that makes the Hamiltonian \f$\Gamma(\mathbf{Q},\mathbf{P})\f$ separable for \f$P\f$ and \f$Q\f$:

(9) \f$ \Gamma(\mathbf{Q},\mathbf{P}) = a(\mathbf{P}) + b(\mathbf{Q}) \f$

Then explicit Leapfrog (SIA) integration method can be used. Preto & Tremaine (1999) suggests to use

(10) \f$ g(\mathbf{Q},\mathbf{P}) = \frac{f(T(\mathbf{P})) - f(-U(\mathbf{Q}))}{T(\mathbf{P}) + U(\mathbf{Q})} \f$

where \f$ T(\mathbf{P}) = T(\mathbf{p}) + Pt \f$ is the extended kinetic energy and \f$ U(\mathbf{Q}) = U(\mathbf{q},t) \f$ is the extended potential energy. 

The Hamiltonian becomes separable:

(11)  \f$ \Gamma = f(T(\mathbf{P})) - f(-U(\mathbf{Q})) \f$

Then the equation of motions are:

(12) \f$ \frac{d \mathbf{q} }{d s} = f'(T(\mathbf{p})+Pt) \frac{\partial T(\mathbf{p})}{\partial {\mathbf{p}}} \f$;
     \f$ \frac{d t }{d s} = f'(T(\mathbf{p})+Pt) \f$;
     \f$ \frac{d \mathbf{p} }{d s} = f'(-U(\mathbf{q},t)) \frac{\partial U(\mathbf{q},t)}{\partial {\mathbf{q}}} \f$;
     \f$ \frac{d Pt}{d s} = f'(-U(\mathbf{q},t)) \frac{\partial U(\mathbf{q},t)}{\partial {\mathbf{t}}} \f$;

where \f$ f'(x) = \frac{d f(x)}{d x} \f$.

Since \f$Pt = -H(t)\f$, \f$H'=H+Pt = T(\mathbf{P}) + U(\mathbf{Q}) = 0 \f$. Thus during integration, \f$T(\mathbf{P}) \approx -U(\mathbf{Q}) \f$. 
This requires \f$ f(T(\mathbf{P})) - f(-U(\mathbf{Q})) \approx 0 \f$. 
With Taylor expansion, we can obtain:

(13) \f$ f(T(\mathbf{P})) = f(-U(\mathbf{Q})) + \left[T(\mathbf{P}) + U(\mathbf{Q})\right] f'(-U(\mathbf{Q})) + O\left[T(\mathbf{P}) + U(\mathbf{Q})\right]^2 \f$

Thus 

(14) \f$ g(\mathbf{Q},\mathbf{P}) \approx f'(-U(\mathbf{Q})) \f$

\subsubsection logH_sec Logarithmic Hamintonian method (LogH)

Mikkola & Tanikawa (1999) and Preto & Tremaine (1999) suggest to use the function \f$ f(x) = \log{x} \f$ (Logarithmic Hamintonian method).
In this case, the time transformation based on (14) is:

(15) \f$ g(\mathbf{Q},\mathbf{P}) \approx \frac{1}{-U(\mathbf{Q})} \f$

Then the equation of motions can be written as:

(16) \f$ \frac{d \mathbf{q} }{d s} = \frac{1}{T(\mathbf{p})+Pt} \frac{\partial T(\mathbf{p})}{\partial {\mathbf{p}}} \f$;
     \f$ \frac{d t }{d s} = \frac{1}{T(\mathbf{p})+Pt} \f$;
     \f$ \frac{d \mathbf{p} }{d s} = \frac{1}{-U(\mathbf{q},t)} \frac{\partial U(\mathbf{q},t)}{\partial {\mathbf{q}}} \f$;
     \f$ \frac{d Pt}{d s} = \frac{1}{-U(\mathbf{q},t)} \frac{\partial U(\mathbf{q},t)}{\partial t} \f$;

For the point mass systems with Newtonian gravity 

(17) \f$ T(\mathbf{p}) = \sum_{i=1}^{n} \frac{\mathbf{p_i}^2}{2m} \f$; \f$ U(\mathbf{q},t) = - \sum_{i<j,i=1,j=1}^{i\rightarrow n,j\rightarrow n} \frac{G m_i m_j}{|\mathbf{q_i}-\mathbf{q_j}|} \f$

where G is gravitational constant and \f$ m_i, m_j \f$ are masses of point-mass particles.

From (17) we see \f$ \frac{d Pt}{d s} = 0 \f$. 
This is only for the isolated system. If the system has external force from perturbers or external potential. The energy of system (\f$-Pt\f$) may not be conserved any more. Thus the energy change should be added into \f$Pt\f$ during the integration.


<A HREF="http://adsabs.harvard.edu/abs/2002CeMDA..84..343M">Mikkola & Aarseth (2002)</A> introduced a modified version of time transformed symplectic method.
Instead of calculating $T(\mathbf{P})$, one can also define a variable,

(18) \f$  u = \int \frac{\partial U(\mathbf{q})}{\partial \mathbf{q}} \cdot \frac{\mathrm{d} \mathbf{q}}{\mathrm{d} t} \f$

Then \f$ f'(u) = f'(T(\mathbf{P})) \f$. The equation of motion has the form:

(19) \f$ \frac{d \mathbf{q} }{d s} = \frac{1}{u} \frac{\partial T(\mathbf{p})}{\partial {\mathbf{p}}} \f$;
     \f$ \frac{d t }{d s} = \frac{1}{u} \f$;
     \f$ \frac{d \mathbf{p} }{d s} = \frac{1}{-U(\mathbf{q})} \frac{\partial U(\mathbf{q})}{\partial {\mathbf{q}}} \f$;
     \f$ \frac{d u}{d s} = \frac{1}{-U(\mathbf{q})} \frac{\partial U(\mathbf{q})}{\partial \mathbf{q}} \cdot \frac{\langle \mathbf{p} \rangle} {m} \f$;

\section code_sec Implementation 

We implememted the SDAR method in this code by using the C++ programming Language. 
This code contains three modules: AR, H4 and COMM.

The AR module provides an c++ integrator class AR::TimeTransformedSymplecticIntegrator.
It is a template class that depends on the class types of particle, interaction, perturbation and extra-information.
AR::TimeTransformedSymplecticManager is the manager class contain the pair interaction class and parameters to control the integration.
One manager can be shared with multiple integrators.

The <A HREF="http://www.sciencedirect.com/science/article/pii/0375960190900923"> Yoshida high-order symplectic method </A> is used together with the SDAR method.
The drift-kick-drift mode must be used as the base of the 2nd order method for a high accuracy.
The reason is described in Preto & Tremaine (1999) and <A HREF="https://github.com/nitadori/Grad4th/blob/master/logH.pdf">Nitadori (2018) </A>.

The Hermite module provides an c++ integrator class H4::HermiteIntegrator that combine a 4th-order Hermite method with block time steps and multiple SDAR integrators.
The former deals with the global particle system and the latter handle the compact subgroups in the system.
The criterion to determine the members in subgroups are based on a distance given by the input parameters.
H4::HermiteManager provides the pair interaction and parameters to control the Hermite integration.

The COMM module provides the basic data type and the tool to construct a Hierarchical (Kepler) binary tree for a group of particles.
It is used to identify the subgroups in the Hermite integrator and to calculate the slow-down factor in the AR integrator.

\section sample_sec Use sample codes

To help the users to understand how to use the library, a few sample codes are provided in the sample directory.
In the three subdirectory: AR, Hermite and Kepler, users can find the instance of code.
For example, AR/ar.cxx is an SDAR integrator, that read a particle set and integrate the system to a given time or a number of steps.
Similarly, Hermite/hermite.cxx is a Hermite integrator.
In Kepler subdirectory, keplerorbit.cxx provides a transformation tool to calculate the Kepler orbital paramters from a particle pair and vice versa.
Keplertree.cxx can construct the Hierarchical Kepler binary tree for a given particle set and vice versa.

Use commander "make" in each subdirectory to create the excutable files and use the sample code.
For example, in AR after make, four ar.** binary files are created.
For any of them (e.g. ar.logh), use 

ar.logh -h

will output the help and describe the input parameters to use the integrator.
Similary, -h options can be used to all other sample codes.
Please refer to the help information for the details of the usage of the sample codes.

\section library_sec For developers

Each c++ class in this library are detailed described in the documantion in order to help the users to develop their own N-body codes.
The developers can use the functions provided in three namespaces: AR, H4 and COMM located in three subdirectories in src.
The sample directory provide the good examples showing how to use the code as a library.

\subsection ar_sec AR module

The AR namespace contains the template classes to construct the SDAR integrator.
The major class that contains the particle data and the integration functions is AR::TimeTransformedSymplecticIntegrator.
It depends on five template types:
- particle type of members
- particle type of the center-of-the-mass of the particle group
- perturber class that handles the external perturbation to the members
- interaction class that defines the pair interaction and the method to calculate slow-down perturbation and timescales
- information class that contains user-defined information, it should inherit from the defaulted AR::Information class.

To create a SDAR integrator, it is necessary to provide these five types.
Thus, the user should first define the particle type.
The example is shown in the class Particle in sample/AR/particle.h.
A few necessary members and member functions should be defined as shown in the sample, please follow the link to check the details.

The particle type of the center-of-the-mass can be different from the type of the member.
One example is shown in the sample/Hermite/hermite.cxx.

If interaction between the members of the particle groups and the external sources is needed.
In the sample of ar.cxx, there is no external perturbation, thus the class Perturber defined in sample/AR/perturber.h is emplty with only a few necessary IO functions

The interaction class is important and necessary to be provided to calculate the pair interaction between the members.
Also the interaction between perturbers and members also also defined in this class.
The example Interaction in sample/AR/interaction.h shows the necessary member functions used in the SDAR method.
For the SDAR method, it is also required to calculate the time tranformation function together with the force and potential calculations.
Besides, if the slow-down method is used, the slow-down inner force and perturbation should also be calculated to calculate the next slow-down factor.

The information class should inherit from AR::Information. 
The latter contain the method to calculate the integration step size based on the Kepler orbital information of the system.
Also it contains the binary tree of the particle members, which records the hierarchical relations between particles and the corresponding Kepler orbital elements.

Once all these 5 types are defined.
Use

    AR::TimeTransformedSymplecticIntegrator<Particle, Particle, Perturber, Interaction, Information<Particle,Particle>> sym_int;
    AR::TimeTransformedSymplecticManager<Interaction> manager;

to create the integrator ARint and the manager.
A few parameters in the manager should be initialized before the integration. Please check the link of AR::TimeTransformedSymplecticManager for details.

There are two way to assign a particle group to the integrator.
use 

    sym_int.particles.setMode([COMM::ListMode::[types]);
    
to set the mode:
- COMM::ListMode::local indicates that the particle data is stored in the local allocated array in particles.
- COMM::ListMode::link indicates that the particle data is not stored in particles but only the address pointing to an existed particle data array is stored.
- COMM::ListMode::copy indicates that the particle data are copied from an existed particle array and also the addresses to original particle are saved in order to writeback the data.

After set mode, the particles can be added based on the mode.
For example

    sym_int.particles.setMode(COMM::ListMode::local);
    sym_int.particles.readMemberAscii(fin);

will read the particle data from a std::fstream file IO fin.

The details of the mode can be found in COMM::List.

The center-of-the-mass can be calculated after all particle data are read.
    
    sym_int.particles.calcCenterOfMass();

Before starting the integration, it is useful to construct the binary tree by using

    sym_int.info.reserveMem(sym_int.particles.getSize());
    sym_int.info.generateBinaryTree(sym_int.particles,manager.interaction.gravitational_constant);

The first line allocate the necessary memory space to store the binary tree.
The second line generate the binary tree. The second argument is the gravitatinal constant, based on the units used in the particle data.

Then the initialization of integration will calculate the initial state of the system (e.g. energy and slow-down factor)

    sym_int.initialIntegration(time_zero.value);
    sym_int.info.calcDsAndStepOption(sym_int.slowdown.getSlowDownFactorOrigin(), manager.step.getOrder(), manager.interaction.gravitational_constant);

The second line estimate the initial step, see AR::Information for details.

Finally, the integration can be excuted by:

    sym_int.integrateToTime(time);

which integrate the system to given time.
If interruption function is triggered, the integration can break in the middle and return the binary tree address of the interrupted pair.
See Interaction class for the interruption function.
Instead of integrating to a given time, it can also only integrate one step by
    
    sym_int.integrateTwoOneStep(sym_int.info.ds, time_table);

or

    sym_int.integrateOneStep(sym_int.info.ds, time_table);

The first is the fast method to integrate two body system. The second is general for multiple systems.

The state of particles can be directly accessed by checking the sym_int.particles or write to std::ostream using the IO functions.
The details of output is defined by users in the writeBinary, readBinary, writeAscii and ReadAscii functions in each class.

\subsection h4_sec H4 module

The way to use Hermite integrator is similar to the SDAR integrator.
The users should provide the additional interaction class which define the interactions between AR members and Hermite members, and the corresponding perturber, information types for Hermite integrator.

To read the particle data, the subgroups that need SDAR methods can be pre-defined in the input data file.
H4::HermiteIntegrator::readGroupConfigureAscii function is used to read the configuration of the subgroups.

The intergration contain four steps (see hermite.cxx sample)

    auto* bin_interrupt = h4_int.integrateGroupsOneStep();
    h4_int.integrateSingleOneStepAct();
    h4_int.adjustGroups(false);
    h4_int.initialIntegration();
    h4_int.sortDtAndSelectActParticle();

The first line integrate all subgroups using SDAR method, if interreupt appears, the address is returned.
In such case, h4_int.integrateGroupsOneStep() should be called again until no new interruption appear.
Then the single particles are integrated by Hermite method.
After one step integration.
H4::HermiteIntegrator::adjustGroups function checks the structure of the system and decide whether subgroups need to form or disrupt.
H4::hermiteIntegrator::initialIntegration renews the system after the adjustment of groups.
The finall step H4::hermiteIntegrator::sortDtAndSelectActParticle sorts the time steps of particles and determine the next active block step lists.

\subsection comm_sec COMM module

In src/Common, a few header files define the common tools used in AR and H4 modules.
- Float.h defines the floating point type, if the QD library is used, high-precision floating point can be switched on.
- list.h define the allocateble array type list, used to manage the particle array.
- binary_tree.h defines the COMM::BinaryTree class that can be used to generate the binary tree of a particle group and COMM::Binary class that can be used to transform particle pair to kepler orbits and vice versa.
- particle_group.h is based on list.h and used as the particle data container for AR and H4 methods.

*/
