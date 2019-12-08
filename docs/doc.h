/*! \mainpage ARC Introduction
  
\section AR_sec Algorithmic Regularization (AR)

The algorithm used in this code is based on the literatures of <A HREF="http://adsabs.harvard.edu/abs/1999MNRAS.310..745M">Mikkola & Tanikawa (1999)</A> and <A HREF="http://adsabs.harvard.edu/abs/1999AJ....118.2532P">Preto & Tremaine (1999)</A>. The development of this code refers to the Chapter 2 (Mikkola) in book <A HREF="http://www.springer.com/us/book/9781402084300">The Cambridge N-body Lectures</A>. Here the basic idea of AR is described.

The numerical simulations of gravitational N-body systems dynamical evolutions are frequently used in astrophysics. 
However, due to the singularity of two-body gravitational potential when these two particles become infinite close lead to the difficulty in the highly accurately integration of two-body bounded system with very high eccentricity. 
To get high accuracy of integration when two particles are very close, the time step for integration should be reduced significantly.
This result in time consuming computation if number of particles is large.
On the other hand, for a long-term integration, the total energy of the systems may be systematiclly drifted due to the numerical accuracy of integrators. 
Thus the sympletic methods are suggested to be used since it can keep the energy conservation for long-term integration.

However, the sympletic methods are difficult to be applied for gravitational systems due to the required time step (integration step) shrinking when two particle get close.
Thus Mikkola & Tanikawa (1999) and Preto & Tremaine (1999) develop the special time transformation method based on the extended phase space Hamiltonian. 
The time \f$t\f$ become a general coordinate in Hamiltonian with corresponding general momentum \f$Pt\f$.
The integration of the equation of motion then depends on the new differential variable \f$ s\f$.
In this case, time and the motion of the system can be integrated with a fixed step size of s, which allow the usage of sympletic methods.

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

Mikkola & Tanikawa (1999) suggests to use the function \f$ f(x) = \log{x} \f$ (Logarithmic Hamintonian method).
In this case, the time transformation based on (14) is:

(15) \f$ g(\mathbf{Q},\mathbf{P}) \approx \frac{1}{-U(\mathbf{Q})} \f$

Then the equation of motions can be written as:

(16) \f$ \frac{d \mathbf{q} }{d s} = \frac{1}{T(\mathbf{p})+Pt} \frac{\partial T(\mathbf{p})}{\partial {\mathbf{p}}} \f$;
     \f$ \frac{d t }{d s} = \frac{1}{T(\mathbf{p})+Pt} \f$;
     \f$ \frac{d \mathbf{p} }{d s} = \frac{1}{-U(\mathbf{q},t)} \frac{\partial U(\mathbf{q},t)}{\partial {\mathbf{q}}} \f$;
     \f$ \frac{d Pt}{d s} = \frac{1}{-U(\mathbf{q},t)} \frac{\partial U(\mathbf{q},t)}{\partial t}} \f$;

For the point mass systems with Newtonian gravity 

(17) \f$ T(\mathbf{p}) = \sum_{i=1}^{n} \frac{\mathbf{p_i}^2}{2m} \f$; \f$ U(\mathbf{q},t) = - \sum_{i<j,i=1,j=1}^{i\rightarrow n,j\rightarrow n} \frac{G m_i m_j}{|\mathbf{q_i}-\mathbf{q_j}|} \f$

where G is gravitational constant and \f$ m_i, m_j \f$ are masses of point-mass particles.

From (17) we see \f$ \frac{d Pt}{d s} = 0 \f$. 
This is only for the isolated system. If the system has external force from perturbers or external potential. The energy of system (\f$-Pt\f$) may not be conserved any more. Thus the energy change should be added into \f$Pt\f$ during the integration.

\subsubsection TTL_sec Time-Transformed Leapfrog (TTL)

The regularization methods where energy explicitly appear in the equation of motions cannot solve the few-body systems with large mass ratio (for example, planetary systems and super massive black hole with surrounding stars), because the energy is dominated by the massive bodies, and this introduce the systematic error during the integration. To solve this kind of issue, <A HREF="http://adsabs.harvard.edu/abs/2002CeMDA..84..343M">Mikkola & Aarseth (2002)</A> developed the so-called Time-Transformed Leapfrog (TTL) method. This method is also based on time transformation. The major difference compared with the LogH method is that the time transformation function also need to be integrated.

The time transformation (10) leads to the equations of motion (12) where time transformation \f$ f'(T(\mathbf{p})+Pt) \f$ and \f$ f'(-U(\mathbf{q},t))\f$ explicitly depend on kinetic energy, binding energy and potential. 
If we want to replace \f$ -U(\mathbf{q},t) \f$ to other quantity \f$ W(\mathbf{q})\f$ (here \f$ W(\mathbf{q})\f$ is positive), considering the requirement \f$ f(T(\mathbf{P})) - f(-U(\mathbf{Q})) \approx 0 \f$, we should also find another quantity \f$ w(\mathbf{p}) \f$ that allow \f$ f(w(\mathbf{p})) - f(W(\mathbf{q})) \approx 0 \f$. and 

(18) \f$ g(\mathbf{Q},\mathbf{P}) = \frac{f(w(\mathbf{p})) - f(-W(\mathbf{q}))}{T(\mathbf{P}) + U(\mathbf{Q})} \approx f'(W(\mathbf{q})) \f$

Instead of finding the \f$ w(\mathbf{p}) \f$ for each kind of \f$ W(\mathbf{q})\f$, Mikkola & Aarseth (2002) suggest to use the differential equation 

(19) \f$  \frac{d W(\mathbf{q})}{d s} = \frac{\partial W(\mathbf{q})}{\partial \mathbf{q}} \cdot \frac{d \mathbf{q}} {d s} \f$

and integrate this equation to approximate \f$ w(\mathbf{p}) = \int \frac{d W(\mathbf{q})}{d s} d s\f$ simultaneously with integration of \f$ \frac{d \mathbf{p} }{d s} \f$.

However 

(20) \f$ \frac{d \mathbf{q}}{d s} = \frac{d \mathbf{q}}{d t} \frac{d t}{d s} = \frac {\mathbf{p}}{m} f'(W(\mathbf{q}))\f$

Thus \f$ \frac{d W(\mathbf{q})}{d s} \f$ explicitly depends on the momemtum. The integration in principle are not separatable. 
To solve this issue, Mikkole & Aarseth (2002) recommend to use averaged momemtums \f$ \langle \mathbf{p} \rangle \f$ (velocities) between previous and current step's during the Leapfrog integration, because the averaged values can represent the momemtums at the D (half) step when \f$\mathbf{q}\f$ is integrated.

Then if we take \f$ f(x) = \log{x}\f$ again, we have the equations of motion like:

(21) \f$ \frac{d \mathbf{q} }{d s} = \frac{1}{w} \frac{\partial T(\mathbf{p})}{\partial {\mathbf{p}}} \f$;
     \f$ \frac{d t }{d s} = \frac{1}{w} \f$;
     \f$ \frac{d \mathbf{p} }{d s} = \frac{1}{W(\mathbf{q})} \frac{\partial U(\mathbf{q},t)}{\partial {\mathbf{q}}} \f$;
     \f$ \frac{d w}{d s} = \frac{1}{W(\mathbf{q})} \frac{\partial W(\mathbf{q})}{\partial \mathbf{q}} \cdot \frac{\langle \mathbf{p} \rangle} {m} \f$;

This solution avoid use the energy (potential) as a time transformation dependence, thus with a suitable choice of \f$ W(\mathbf{q}) \f$, the high mass ratio systems can be integrated with high accuracy.

\section code_sec Implementation of ARC

We implememted AR method together with Chain (discussed below) for few-body systems by using C++ programming Language. 
The idea is make the integrator a C++ class thus can be easily used as a module for other codes. 
In this section we describe the details of the implementation.

\subsection chain_sec Particle Chain

If the bounded few-body systems are inside a big cluster enviroment, the average distance between these particles can be much smaller than the scale of cluster. 
Thus the round off error can be large if the positions of these particles are in the cluster center-of-mass frame.
To avoid this issue, <A HREF="http://adsabs.harvard.edu/abs/1993CeMDA..57..439M">Mikkola & Aarseth (1993)</A> suggested to use Chain method.

The idea is to connect all particles in one chain and using relative position and velocity for integration.
Firstly, one particle is selected as a starting point of the chain, then the nearest particle is selected as the next chain member, the relative position \f$ X \f$ and velocity \f$ V \f$ between these neighbors are calculated and stored.
After that, we found the nearest particle to this second member from the remaining particles and calculate relative positions and velocites and do this iterately until all particles are connected in this chain.
The relative positions and velocites can be described by absolute positions and velocities in a ordered chain as:

(22) \f$ \mathbf{X}_i = \mathbf{q}_{i+1} - \mathbf{q}_i \f$; \f$ \mathbf{V}_i = \mathbf{v}_{i+1} - \mathbf{v}_i \f$

The integration is done with these relative quantities to reduce round off error. The equations of motion can be written as 

(23) \f$ \frac{d \mathbf{X}_i}{d t} = \mathbf{V}_i \f$; \f$ \frac{d \mathbf{V}_i}{d t} = \mathbf{A}_{i+1} - \mathbf{A}_i \f$

where \f$ \mathbf{A}_i \f$ is the acceleration of particle \f$ i\f$.

When the particles are moved, the nearest neighbor of each particle may become different, thus the update of chain order should be performed with a suitable time interval. 

\subsection leap_sec Leapfrog Integrator

By combining the AR algorithm and Chain scheme, we can construct a Leapfrog integrator of equations of motion for $N$-body systems like:
- D mode:

(24)      \f$ \Delta t = \Delta s / (\alpha (T(\mathbf{p}) + Pt) + \beta w + \gamma) \f$;
          \f$ t += \Delta t \f$;
          \f$ \mathbf{X}_i += \Delta t \mathbf{V}_i \f$ 

- K mode:

(25)      \f$ \delta t = \Delta s / (\alpha U(\mathbf{q},t) + \beta W(\mathbf{q}) + \gamma) \f$;
          \f$ \mathbf{V}_i += \delta t (\mathbf{A}_{i+1} - \mathbf{A}_{i}) \f$;
          \f$ Pt += \delta t \sum_i (-m_i \langle \mathbf{v}_i \rangle \cdot f_{ext,i}) \f$;
          \f$ w += \delta t \sum_i \frac{\partial W}{\partial \mathbf{q}_i} \cdot \langle \mathbf{v}_i \rangle \f$ 

where \f$ f_{ext,i} \f$ is the external force from outside the system (e.g., perturber force or tidal force) of each particle \f$ i\f$, and \f$ \langle \mathbf{v}_i \rangle\f$ is obtained by averaging the velocities of the initial and the final \f$ \mathbf{v}_i \f$ of this K mode step. \f$ \alpha, \beta, \gamma \f$ are the coefficients representing the weights of the LogH, TTL and non-time-transformation modes separately. For example, if \f$ \alpha=0\f$, then no LogH will be performed, and if \f$ \alpha =1, \beta=0, \gamma=0 \f$ it is LogH ARC.

The initial value of \f$ Pt \f$ should be the initial binding energy of the system \f$ U(\mathbf{q},t) - T(\mathbf{p}) \f$. 
If the system is isolated, \f$ Pt \f$ is constant.
The initial value of \f$ w\f$ is set to initial \f$ W(\mathbf{q}) \f$.

The Leapfrog step start with half-step D and then loop full-step K-D-K and stop with half-step D:

(26) \f$ D(\Delta s/2)K(\Delta s)D(\Delta s)....K(\Delta s)D(\Delta s/2) \f$

This provide a second order integrator of ARC. Trying this integrator for a two-body bounded system can result in an energy and eccentricity conserved kepler orbit. Only the time phase can have cumulative error after long-term integration.

\subsection extrapolation_sec Extrapolation Integrator

The Leapfrog integrator only has second order accuracy, which is not enough for many applications. One can reduce the step size of integration to obtain higher accuracy. 
However, as energy is always conserved for two-body motions, we don't have good checker to indicate whether the integration is accurate enough. 
A better and more efficient way is to extrapolate the integration results to infinite small step \f$ \Delta s\approx 0\f$, thus the high accuracy result can be obtained.
The idea of extrapolation integration is well summarized in <A HREF="http://link.springer.com/book/10.1007%2F978-0-387-21738-3"> Stoer & Bulirsch</A>. 
Here the basic algorithm is shown.

First, if we integrate the equations of motion with Leapfrog integrator by step \f$ \Delta s\f$. 
we get the first result with a certain accuracy. Now we keep the total step constant but divide the integration into several sub-steps with equal sizes by \f$ n \f$, we can obtain higher accuracy of the integration. When we use a sequence of dividers \f$ (n_1, n_2, n_3 ...)\f$ (\f$ n_{i+1}>n_i\f$) and do the integration with each of them, we can obtain a series of results with increasing accuracy. Then we can extrapolate these results to obtain the value of \f$ \Delta s/n_{\infty}=0 \f$.

There are two major methods of extrapolation: polynomial and rational.
Both methods can be described as recursive functions:

- Polynomial: 

(27) \f$ T_{i,k} = T_{i,k-1} + \frac{T_{i,k-1} - T_{i-1,k-1}}{( h_{i-k} / {h_i} )^2 -1} \f$, \f$ 1 \le k \le i \le m \f$ 

- Rational:

(28) \f$ T_{i,k} = T_{i,k-1} + \frac{T_{i,k-1} - T_{i-1,k-1}}{\left[ \frac{h_{i-k}}{h_i} \right]^2 \left[ 1 - \frac{T_{i,k-1} - T_{i-1,k-1}}{T_{i,k-1}- T_{i-1,k-2}} \right]-1} \f$, \f$ 1 \le k \le i \le m \f$ 

Here \f$ i\f$ indicate the integration with sub-step size \f$ h_i = s/n_i\f$, and \f$ k \f$ indicate the extrapolation order.
The \f$ T_{i,0} \f$ are results of Leapfrog integrations, and for each order \f$ i\f$, the \f$ T_{i,i} \f$ is final extrapolation result we want. 
The \f$ T_{i,i} \f$ can be obtained by calculating \f$ T_{i,k} \f$ from \f$ k=1 \f$ to \f$ k=i \f$ using the recursive functions.

One benefit of these recursive functions is that a higher order extrapolation \f$ T_{i+1,i+1} \f$ can be established based on current existing \f$ T_{i,k}, k=0\sim i \f$  with a new higher order integration result \f$ T_{i+1,0} \f$. 
Then it is easy to estimate the error by comparing \f$ T_{i+1,i+1} \f$ and \f$ T_{i,i} \f$ to determine whether another higher order result is necessary. 
For example, in ARC integration, we can check the time or position phase error and energy error to determine how many orders we need to integrate and extrapolate due to the accuracy requirment.

The sequences of dividers \f$ n_i \f$ have several choices for different applications:
- <A HREF="https://en.wikipedia.org/wiki/Romberg's_method">Romberg</A>: (1, 2, 4, 8 ...)
- <A HREF="http://link.springer.com/article/10.1007%2FBF02165234">Bulirsch & Stoer</A> (BS): (1, 2, 3, 4, 6, 8 ...)
- <A HREF="http://link.springer.com/article/10.1007%2FBF01385634">Hairer</A> (4k): (2, 6, 10, 14 ...)
- Harmonic:  (1, 2, 3, 4 ...)

Different seuqnces and recursive functions can be combined together for extrapolation integration. 
We implement all sequences shown above. Later we discuss the special application of some sequences.

\subsection dense_sec Dense Output for Time Synchronization

Although the ARC can make the integration of $N$-body systems accurately, the side-effect of time transformation is that the physical time become unpredictable.
With the Leapfrog integrator, we cannot know what will be the final physical time before one integration step finish.
This result in difficulty if we want to use the ARC together with a $N$-body code to simulate a particle cluster including dense sub-systems.
The integration of the motions of particles surrounding this sub-system need to obtain the acceleration from this sub-system at a certain physical time, but with ARC the integration of this sub-system cannot exactly reach the required time.
Especially with extrapolation method, the large integration step is used frequently, thus the physical time error can be significant.

To solve this issue, we apply the dense output of extrapolation method introduced by <A HREF="http://link.springer.com/article/10.1007%2FBF01385634">Hairer & Ostermann (1990)</A>. 
The idea of this scheme is using interpolation to obtain the integrated variable at any sub-step inside an extrapolation integration step.
The interpolation should have the similar order of accuracy as the extrapolation and the internal integration results during extrapolation should be used for interpolation to save computational effort.

The physical time \f$ t\f$ as a function of integration step variable \f$ s \f$ then can be interpolated as \f$ T(s) \f$. 
If the required ending physical time \f$ t_{off} \f$ is inside one integration step, we can solve the equation \f$ T(s)=t_{off} \f$ to obtain the correct step size \f$\Delta s_{off} \f$ to reach the exact \f$ t_{off} \f$.
Then by redoing this integration step with \f$\Delta s_{off}\f$, we can get correct results.

One can also try to do dense ouput for all variables (\f$t\f$, \f$Pt\f$, \f$w\f$, \f$\mathbf{q}\f$, \f$\mathbf{p}\f$), thus the results at correct physical time can be directly calculated instead of redoing the integration.
However, as the computation of dense output is quite heavy (many extrapolation is needed; see below), redoing the integration can be cheaper if particle number is not large (\f$<=4\f$).


Hairer & Ostermann (1990) introduced two dense output methods.
One is for explicit Euler integrator using Harmonic sequences and another is for Gragg-Bulirsch-Stoer (GBS) method with 4k sequences (shown above).
Here the brief algorithms are shown without mathematical proof.

\subsubsection euler_dense_sec Dense Output for explicit Euler

If the integrated variable is \f$ y \f$ and its first derivate (acceleration) is \f$ f \f$ which can be calculated directly, we can use explicit Euler together with Harmonic sequence \f$ n_i = i\f$ for extrapolation.
Then during each integration step, we have the initial \f$ y_i(0) \f$ and the final \f$ y_i(\Delta s) \f$.
In addition, \f$ f_i(\Delta s* k/n_i) (k=0,n_i)\f$ are also calculated.
Thus we can obtain the high order derivates of \f$ f_i \f$ at the left and right edges using forward and backward differences:

(29) \f$ f_i^{(k)}(0) = \left[\frac{n_i}{\Delta s}\right]^k \sum_{j=0}^k (-1)^j B_j^k f( \Delta s*\frac{k-j}{n_i}) \f$;  \f$ f_i^{(k)}(\Delta s) = \left[\frac{n_i}{\Delta s}\right]^k \sum_{j=0}^k (-1)^j B_j^k f(\Delta s*(1-\frac{j}{n_i})) \f$; \f$(k=1, n_i) \f$

where \f$ B_j^k = \frac{j!}{k!(j-k)!}\f$ is the binomial sequence.

Then if the last sequence index used in extrapolation is \f$i=\kappa\f$, the maximum order of derivate is \f$ n_\kappa\f$.
Besides, for each order of derivate \f$ f_i^{(k)}(0) \f$ and \f$ f_i^{(k)}(\Delta s)\f$ (\f$ k=1,n_\kappa\f$), we also have the values of different order of accuracy corresponding to difference step size \f$(\Delta s/n_i)\f$ (\f$ i=k,n_\kappa\f$). 
Thus the extrapolation can be done with these different order of accuracy (the same way as the extrapolation of \f$y(\Delta s)\f$) to get high accurate derivates \f$f^{(k)}(0)\f$ and \f$f^{(k)}(\Delta s)\f$.

Since now the \f$ y(0)\f$, \f$y(\Delta s)\f$, \f$f^{(k)}(0)\f$ and \f$f^{(k)}(\Delta s)\f$ are avaiable, then Hermite interpolation can be used to get the interpolation polynomial function \f$ Y(x) \f$ and

(30) \f$ Y(x) - y(x) = O(\Delta s^{n_\kappa+1}) \f$

where \f$ n_\kappa = \kappa \f$ in the case of Harmonic sequence.

\subsubsection gbs_dense_sec Dense Output for Gragg-Bulirsch-Stoer

Similar as the dense output method described above, for mordified middle point integrator used in GBS method, we can construct the interpolation using high order derivates of \f$ f \f$ at the middle position (\f$\Delta s/2\f$) instead of edges.
However, differing from the edge differences, the middle difference is sensitive to the data point number.
If \f$ n_i\f$ is even, to obtain the derivate order with odd \f$ k \f$ (which means \f$k+1\f$ points are needed), we have to use values every two sub-steps.

For example, when \f$ n_i = 6 \f$, there are 6 sub-steps and 7 points (\f$ \Delta s*j/n_i \f$ with \f$j=0,6\f$).
If \f$ k = 3 \f$, 4 points are needed to obtain the derivate \f$ f_i^{(3)}(\Delta s/2) \f$. 
Since we need the derivate at \f$ \Delta s/2 \f$, only \f$ f \f$ at \f$ j = 0,2,4,6 \f$ can be used.
If \f$ k = 4 \f$, values at \f$ j= 1,2,3,4,5\f$ are OK.
But the difference step sizes in this case are different for odd and even \f$ k\f$.

To keep accuracy order consistent, we only allow every two points to be used for both odd and even order of derivates.
The formular then should be

(31) \f$ f_i^{(k)}(\Delta s/2) = \left[ \frac{2n_i}{\Delta s} \right]^k \sum_{j=0}^k (-1)^j B_j^k f(\Delta s*(\frac{1}{2}+\frac{z_j-2j}{n_i})) \f$; \f$ k=1,2i-1 \f$
- if \f$k\f$ is odd, \f$ z_j = k+1 \f$
- if \f$k\f$ is even, \f$ z_j = k \f$

together with the 4k sequence \f$ n_i =(2, 6, 10, 14 ...) \f$.

Then again the extrapolation of the derivates and also the middle point integrated variable \f$ y(\Delta s/2) \f$ can be done and \f$ y(0) \f$, \f$ y(\Delta s) \f$, \f$ y(\Delta s/2) \f$ and derivates \f$ f^{(k)}(\Delta s/2) \f$ (\f$ k =1,2\kappa-1 \f$) are avaiable for Hermite interpolation.
This method can provide the interpolation polynomial function with accuracy

(32) \f$ Y(x) - y(x) = O(\Delta s^{2\kappa-1}) \f$

\subsubsection gbs_dense_instab Instability Issue

Sometimes, when the integrated step size \f$ \Delta s \f$ is large, the interpolation may fails near the edges of the interval. 
Especially when the maximum derivative order is high (>10), the instability may happen and create sharp peaks around the starting and ending points.
To avoid this, the integrated step size should be controlled.
Since there is no theoretical formula to define what is a good step size, the user should try to adjust the step size depending on the problem to solve.
One method that can help to reduce the error is to include edge derivatives to improve the interpolation accuracy around the edges. 
But this require more computational effort and cannot remove this issue completely.

If the variables depending on \f$ s \f$ is monotonic, to detect whether an instability happen, a monotonic test can be performed when using the interpolating polynomial function. 
Since physical time monotonically depends on \f$ s \f$, this method is very useful to avoid serious errors.

\subsection step_sec Integration Step Control

If we use the automatical accuracy order in extrapolation integration (the maximum sequence index \f$\kappa \f$ is determined by the error criterion), the step size \f$\Delta s\f$ can be constant with a suitable initial value.
On the other hand, \f$ \Delta s\f$ can be also adjusted based on integration error to approach better performance.

\subsubsection step_error Step Estimation Based On Extrapolation Error

The integration error during extrapolation at sequence index \f$ i\f$ can be estimated as

(33) \f$ err_i = \frac{2|T_{i,i-1} - T_{i,i}|}{\sqrt{T_{i,i-1}^2 + T_{i,i}^2}}\f$ 

If we want the expected error appear at sequence index \f$ i\f$ after the next integration step, the step modification factor can be estimated as:

(34) \f$ \frac{\Delta s_{new,i}}{\Delta s} \approx \left(\frac{exp}{err_{i}}\right)^{1/(2i-1)} \f$

with the assumption \f$ err_i \propto (\Delta s)^{2i-1} \f$. 
To determine which \f$ i\f$ is best for performance, the computational effort 

(35) \f$ C_i = \frac{\sum_{k=0}^i n_i}{\Delta s_{new,i}} \f$

is calculated for each \f$ i \f$, then we choose index \f$ i=k \f$ which corresponds to the mimimum \f$ C_i \f$.
The next step is \f$ \Delta s_{new,i} \f$.

This method should work with fixed accuracy order (\f$ \kappa \f$ is constant).
But in some critial situation (close encounter), the step change may not be reduced enough to obtain accurate results. 
This need to be treated carefully during the simulations.
On the other hand, in some special cases that the integration cannot reach the energy criterion when reducing the step size, this algorithm will lead to an continuing decreasing of step sizes and significantly influence the performance.
Thus the user should be very careful to use this method.

\subsubsection step_kepler Step Estimation Based On Kepler Period

If we assume every neighbor members in the chain are two-body systems, we can calculate the Kepler period of each pair.
Then the mimimum period can be used to estimate the next step size (a few steps should be carried out for one period).
This may fail if the systems are chaotic or suffering hyperbolic encounters.

In the case of hyperbolic encounters, the free-fall time scale 

(36) \f$ t_f = \frac{\pi}{2} \sqrt{\frac{|\mathbf{X}|^3}{2 m}} \f$

can be used to estimate the next step.
However, the instability of interpolation during dense output can easily happen with this time estimation.
Users should be very careful when treat hyperbolic encounters.

\subsection perf_sec Performance Analysis

Here the performance analysis of the code is provided.
For one step of Leapfrog integration, we need two half-step integration of \f$ \mathbf{q} \f$ and \f$ t \f$, one full-step integration of \f$ \mathbf{p} \f$, \f$ Pt \f$ and \f$ w \f$.
Before \f$ \mathbf{p} \f$ is integrated, the acceleration \f$ \mathbf{A} \f$ is calculated.
If the particle number is \f$ N \f$, the computational cost is:

(37) \f$ C_{LF} = C_{A,P}*N^2 + C_{Pt}*N + C_{w}*N +  2C_{p}*N + 2*C_{t} + C_{T}*N \f$

where \f$ C_* \f$ correspond to the number of operations of different parts. If there is no perturbation and external force, \f$ Pt \f$ is constant and \f$ C_{Pt} = 0\f$. If TTL method is switched off, \f$ C_{w} = 0\f$.

During the extrapolation integration, the Leapfrog integration is performed many times. After integration finished at each sequence \f$ n_i \f$, the extrapolation is performed.
Thus the total cost is:

(38) \f$ C_{EINT} = \sum_{i=1}^\kappa n_i*(C_{LF} + C_{EX}*(6N+3))\f$

where \f$ C_{EX} \f$ is the number of operations of (polynomial or rational) extrapolation function. The \f$ (6N+3) \f$ includes the variables of \f$ t \f$, \f$ Pt \f$, \f$ w \f$, \f$ \mathbf{q} \f$ and \f$ \mathbf{p} \f$.

For the dense output, the high order derivates of \f$ dt/ds \f$ and their extrapolation are calculated.
The cost is:

(39) \f$ C_{DEN} = \sum_{i=1}^\kappa \sum_{j=1}^{2i-1} [2j*C_{DIFF} +  C_{EX}] \f$

where \f$ C_{DIFF} \f$ is the number of operations for adding one \f$ f(x) \f$ value during the computation of difference (31). 
For the two dense output methods discussed Section \ref dense_sec, the cost formula is similar (but the \f$ \kappa \f$ can be significant difference in practice).

As we discussed in Section \ref dense_sec, we can do interpolation for all variables and the cost of dense output is \f$ C_{DEN}*(6N+3) \f$. 
Then the cost of dense output over extrapolation integration is

(40) \f$ \frac{C_{DEN}}{C_{EINT}} \approx \frac{O(\kappa^3)}{O(\langle n_i\rangle*N)} \f$

where \f$\langle n_i\rangle\f$ is the average \f$ n_i \f$ from \f$ i=1,\kappa\f$.
In the case of 4k sequence, \f$\langle n_i\rangle  \propto \kappa\f$.
The value of \f$\kappa\f$ depends on the computational error criterion and the integration step size \f$ \Delta s\f$. 
Usually \f$ \kappa>4 \f$, thus if \f$ N \f$ is not large (\f$ N < 5 \f$), the full dense output with all variables can be more computational expensive.

*/
