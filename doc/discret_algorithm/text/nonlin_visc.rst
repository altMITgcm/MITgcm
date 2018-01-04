Nonlinear Viscosities for Large Eddy Simulation
===============================================

In Large Eddy Simulations (LES), a turbulent closure needs to be
provided that accounts for the effects of subgridscale motions on the
large scale. With sufficiently powerful computers, we could resolve the
entire flow down to the molecular viscosity scales
(:math:`L_{\nu}\approx 1 \rm cm`). Current computation allows perhaps
four decades to be resolved, so the largest problem computationally
feasible would be about 10m. Most oceanographic problems are much larger
in scale, so some form of LES is required, where only the largest scales
of motion are resolved, and the subgridscale’s effects on the
large-scale are parameterized.

To formalize this process, we can introduce a filter over the
subgridscale L: :math:`u_\alpha\rightarrow \BFKav u_\alpha` and
:math:`L:
b\rightarrow \BFKav b`. This filter has some intrinsic length and time
scales, and we assume that the flow at that scale can be characterized
with a single velocity scale (:math:`V`) and vertical buoyancy gradient
(:math:`N^2`). The filtered equations of motion in a local Mercator
projection about the gridpoint in question (see Appendix for notation
and details of approximation) are:

.. math::

   \begin{aligned}
   \BFKaDt \BFKatu- \frac{\BFKatv
     \sin\theta}{\BFKRo\sin\theta_0}+\frac{\BFKMr}{\BFKRo}\BFKpd{x}{\BFKav\pi}
   & = & -\left({\BFKav{\BFKDt \BFKtu}}-{\BFKaDt \BFKatu}\right)
   +\frac{\nabla^2{\BFKatu}}{\BFKRe}\label{eq:mercat}\\
   \BFKaDt\BFKatv+\frac{\BFKatu\sin\theta}{\BFKRo\sin\theta_0}
   +\frac{\BFKMr}{\BFKRo}\BFKpd{y}{\BFKav\pi}
   & = & -\left({\BFKav{\BFKDt \BFKtv}}-{\BFKaDt \BFKatv}\right)
   +\frac{\nabla^2{\BFKatv}}{\BFKRe}\nonumber\\
   \BFKaDt {\BFKav w} +\frac{\BFKpd{z}{\BFKav\pi}-\BFKav b}{\BFKFr^2\lambda^2}
   & = & -\left(\BFKav{\BFKDt w}-\BFKaDt {\BFKav{w}}\right)
   +\frac{\nabla^2\BFKav w}{\BFKRe}\nonumber\\
   \BFKaDt{\ \BFKav b}+\BFKav w & = &
    -\left(\BFKav{\BFKDt{b}}-\BFKaDt{\ \BFKav b} \right)
   +\frac{\nabla^2 \BFKav b}{\Pr\BFKRe}\nonumber \\
   \mu^2\left(\BFKpd x\BFKatu  + \BFKpd y\BFKatv \right)+\BFKpd z {\BFKav w} 
   & = & 0\label{eq:cont}\end{aligned}

 Tildes denote multiplication by :math:`\cos\theta/\cos\theta_0` to
account for converging meridians.

The ocean is usually turbulent, and an operational definition of
turbulence is that the terms in parentheses (the ’eddy’ terms) on the
right of ([eq:mercat]) are of comparable magnitude to the terms on the
left-hand side. The terms proportional to the inverse of , instead, are
many orders of magnitude smaller than all of the other terms in
virtually every oceanic application.

Eddy Viscosity
--------------

A turbulent closure provides an approximation to the ’eddy’ terms on the
right of the preceding equations. The simplest form of LES is just to
increase the viscosity and diffusivity until the viscous and diffusive
scales are resolved. That is, we approximate:

.. math::

   \begin{aligned}
   \left({\BFKav{\BFKDt \BFKtu}}-{\BFKaDt \BFKatu}\right)
   \approx\frac{\nabla^2_h{\BFKatu}}{\BFKRe_h}
   +\frac{\BFKpds{z}{\BFKatu}}{\BFKRe_v}\label{eq:eddyvisc}, & &
   \left({\BFKav{\BFKDt \BFKtv}}-{\BFKaDt \BFKatv}\right)
   \approx\frac{\nabla^2_h{\BFKatv}}{\BFKRe_h}
   +\frac{\BFKpds{z}{\BFKatv}}{\BFKRe_v}\nonumber\\
   \left(\BFKav{\BFKDt w}-\BFKaDt {\BFKav{w}}\right)
   \approx\frac{\nabla^2_h\BFKav w}{\BFKRe_h}
   +\frac{\BFKpds{z}{\BFKav w}}{\BFKRe_v}\nonumber, & &
   \left(\BFKav{\BFKDt{b}}-\BFKaDt{\ \BFKav b} \right)
   \approx\frac{\nabla^2_h \BFKav b}{\Pr\BFKRe_h}
   +\frac{\BFKpds{z} {\BFKav b}}{\Pr\BFKRe_v}\nonumber\end{aligned}

Reynolds-Number Limited Eddy Viscosity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One way of ensuring that the gridscale is sufficiently viscous (*ie.*
resolved) is to choose the eddy viscosity :math:`A_h` so that the
gridscale horizontal Reynolds number based on this eddy viscosity,
:math:`\BFKRe_h`, to is O(1). That is, if the gridscale is to be
viscous, then the viscosity should be chosen to make the viscous terms
as large as the advective ones. Bryan *et al*
:raw-latex:`\cite{Bryanetal75}` notes that a computational mode is
squelched by using :math:`\BFKRe_h<`\ 2.

MITgcm users can select horizontal eddy viscosities based on
:math:`\BFKRe_h` using two methods. 1) The user may estimate the
velocity scale expected from the calculation and grid spacing and set
the viscAh to satisfy :math:`\BFKRe_h<2`. 2) The user may use
viscAhReMax, which ensures that the viscosity is always chosen so that
:math:`\BFKRe_h<{\sf viscAhReMax}`. This last option should be used with
caution, however, since it effectively implies that viscous terms are
fixed in magnitude relative to advective terms. While it may be a useful
method for specifying a minimum viscosity with little effort, tests
:raw-latex:`\cite{Bryanetal75}` have shown that setting viscAhReMax=2
often tends to increase the viscosity substantially over other more
’physical’ parameterizations below, especially in regions where
gradients of velocity are small (and thus turbulence may be weak), so
perhaps a more liberal value should be used, *eg.* viscAhReMax=10.

While it is certainly necessary that viscosity be active at the
gridscale, the wavelength where dissipation of energy or enstrophy
occurs is not necessarily :math:`L=A_h/U`. In fact, it is by ensuring
that the either the dissipation of energy in a 3-d turbulent cascade
(Smagorinsky) or dissipation of enstrophy in a 2-d turbulent cascade
(Leith) is resolved that these parameterizations derive their physical
meaning.

Vertical Eddy Viscosities
~~~~~~~~~~~~~~~~~~~~~~~~~

Vertical eddy viscosities are often chosen in a more subjective way, as
model stability is not usually as sensitive to vertical viscosity.
Usually the ’observed’ value from finescale measurements, etc., is used
(*eg.* viscAr\ :math:`\approx1\times10^{-4} m^2/s`). However,
Smagorinsky :raw-latex:`\cite{Smagorinsky93}` notes that the Smagorinsky
parameterization of isotropic turbulence implies a value of the vertical
viscosity as well as the horizontal viscosity (see below).

Smagorinsky Viscosity
~~~~~~~~~~~~~~~~~~~~~

Some :raw-latex:`\cite{sm63,Smagorinsky93}` suggest choosing a viscosity
that *depends on the resolved motions*. Thus, the overall viscous
operator has a nonlinear dependence on velocity. Smagorinsky chose his
form of viscosity by considering Kolmogorov’s ideas about the energy
spectrum of 3-d isotropic turbulence.

Kolmogorov suppposed that is that energy is injected into the flow at
large scales (small :math:`k`) and is ’cascaded’ or transferred
conservatively by nonlinear processes to smaller and smaller scales
until it is dissipated near the viscous scale. By setting the energy
flux through a particular wavenumber :math:`k`, :math:`\epsilon`, to be
a constant in :math:`k`, there is only one combination of viscosity and
energy flux that has the units of length, the Kolmogorov wavelength. It
is :math:`L_\epsilon(\nu)\propto\pi\epsilon^{-1/4}\nu^{3/4}` (the
:math:`\pi` stems from conversion from wavenumber to wavelength). To
ensure that this viscous scale is resolved in a numerical model, the
gridscale should be decreased until :math:`L_\epsilon(\nu)>L` (so-called
Direct Numerical Simulation, or DNS). Alternatively, an eddy viscosity
can be used and the corresponding Kolmogorov length can be made larger
than the gridscale,
:math:`L_\epsilon(A_h)\propto\pi\epsilon^{-1/4}A_h^{3/4}` (for Large
Eddy Simulation or LES).

There are two methods of ensuring that the Kolmogorov length is resolved
in MITgcm. 1) The user can estimate the flux of energy through spectral
space for a given simulation and adjust grid spacing or viscAh to ensure
that :math:`L_\epsilon(A_h)>L`. 2) The user may use the approach of
Smagorinsky with viscC2Smag, which estimates the energy flux at every
grid point, and adjusts the viscosity accordingly.

Smagorinsky formed the energy equation from the momentum equations by
dotting them with velocity. There are some complications when using the
hydrostatic approximation as described by Smagorinsky
:raw-latex:`\cite{Smagorinsky93}`. The positive definite energy
dissipation by horizontal viscosity in a hydrostatic flow is
:math:`\nu D^2`, where D is the deformation rate at the viscous scale.
According to Kolmogorov’s theory, this should be a good approximation to
the energy flux at any wavenumber :math:`\epsilon\approx\nu D^2`.
Kolmogorov and Smagorinsky noted that using an eddy viscosity that
exceeds the molecular value :math:`\nu` should ensure that the energy
flux through viscous scale set by the eddy viscosity is the same as it
would have been had we resolved all the way to the true viscous scale.
That is, :math:`\epsilon\approx
A_{hSmag} \BFKav D^2`. If we use this approximation to estimate the
Kolmogorov viscous length, then

.. math::

   L_\epsilon(A_{hSmag})\propto\pi\epsilon^{-1/4}A_{hSmag}^{3/4}\approx\pi(A_{hSmag}
   \BFKav D^2)^{-1/4}A_{hSmag}^{3/4} = \pi A_{hSmag}^{1/2}\BFKav D^{-1/2}

 To make :math:`L_\epsilon(A_{hSmag})` scale with the gridscale, then

.. math:: A_{hSmag} = \left(\frac{{\sf viscC2Smag}}{\pi}\right)^2L^2|\BFKav D|

 Where the deformation rate appropriate for hydrostatic flows with
shallow-water scaling is

.. math::

   |\BFKav D|=\sqrt{\left(\BFKpd{x}{\BFKav \BFKtu}-\BFKpd{y}{\BFKav \BFKtv}\right)^2
   +\left(\BFKpd{y}{\BFKav \BFKtu}+\BFKpd{x}{\BFKav \BFKtv}\right)^2}

 The coefficient viscC2Smag is what an MITgcm user sets, and it replaces
the proportionality in the Kolmogorov length with an equality. Others
:raw-latex:`\cite{griffies:00}` suggest values of viscC2Smag from 2.2 to
4 for oceanic problems. Smagorinsky :raw-latex:`\cite{Smagorinsky93}`
shows that values from 0.2 to 0.9 have been used in atmospheric
modeling.

Smagorinsky :raw-latex:`\cite{Smagorinsky93}` shows that a corresponding
vertical viscosity should be used:

.. math::

   A_{vSmag}=\left(\frac{{\sf viscC2Smag}}{\pi}\right)^2H^2
   \sqrt{\left(\BFKpd{z}{\BFKav \BFKtu}\right)^2
   +\left(\BFKpd{z}{\BFKav \BFKtv}\right)^2}

 This vertical viscosity is currently not implemented in MITgcm
(although it may be soon).

Leith Viscosity
~~~~~~~~~~~~~~~

Leith :raw-latex:`\cite{Leith68,Leith96}` notes that 2-d turbulence is
quite different from 3-d. In two-dimensional turbulence, energy cascades
to larger scales, so there is no concern about resolving the scales of
energy dissipation. Instead, another quantity, enstrophy, (which is the
vertical component of vorticity squared) is conserved in 2-d turbulence,
and it cascades to smaller scales where it is dissipated.

Following a similar argument to that above about energy flux, the
enstrophy flux is estimated to be equal to the positive-definite
gridscale dissipation rate of enstrophy :math:`\eta\approx A_{hLeith}
|\nabla\BFKav \omega_3|^2`. By dimensional analysis, the
enstrophy-dissipation scale is :math:`L_\eta(A_{hLeith})\propto\pi
A_{hLeith}^{1/2}\eta^{-1/6}`. Thus, the Leith-estimated length scale of
enstrophy-dissipation and the resulting eddy viscosity are

.. math::

   \begin{aligned}
   L_\eta(A_{hLeith})\propto\pi A_{hLeith}^{1/2}\eta^{-1/6}
   & = & \pi A_{hLeith}^{1/3}|\nabla \BFKav \omega_3|^{-1/3} \\
   A_{hLeith} & = & 
   \left(\frac{{\sf viscC2Leith}}{\pi}\right)^3L^3|\nabla \BFKav\omega_3| \\
   |\nabla\omega_3| & \equiv & 
   \sqrt{\left[\BFKpd{x}{\ }
       \left(\BFKpd{x}{\BFKav \BFKtv}-\BFKpd{y}{\BFKav
           \BFKtu}\right)\right]^2
     +\left[\BFKpd{y}{\ }\left(\BFKpd{x}{\BFKav \BFKtv}
         -\BFKpd{y}{\BFKav \BFKtu}\right)\right]^2}\end{aligned}

Modified Leith Viscosity
~~~~~~~~~~~~~~~~~~~~~~~~

The argument above for the Leith viscosity parameterization uses
concepts from purely 2-dimensional turbulence, where the horizontal flow
field is assumed to be divergenceless. However, oceanic flows are only
quasi-two dimensional. While the barotropic flow, or the flow within
isopycnal layers may behave nearly as two-dimensional turbulence, there
is a possibility that these flows will be divergent. In a
high-resolution numerical model, these flows may be substantially
divergent near the grid scale, and in fact, numerical instabilities
exist which are only horizontally divergent and have little vertical
vorticity. This causes a difficulty with the Leith viscosity, which can
only responds to buildup of vorticity at the grid scale.

MITgcm offers two options for dealing with this problem. 1) The
Smagorinsky viscosity can be used instead of Leith, or in conjunction
with Leith–a purely divergent flow does cause an increase in Smagorinsky
viscosity. 2) The viscC2LeithD parameter can be set. This is a damping
specifically targeting purely divergent instabilities near the
gridscale. The combined viscosity has the form:

.. math::

   \begin{aligned}
   A_{hLeith} & = & 
   L^3\sqrt{\left(\frac{{\sf viscC2Leith}}{\pi}\right)^6
     |\nabla \BFKav \omega_3|^2
     +\left(\frac{{\sf viscC2LeithD}}{\pi}\right)^6
     |\nabla \nabla\cdot \BFKav {\tilde u}_h|^2} \\
   |\nabla \nabla\cdot \BFKav {\tilde u}_h| & \equiv & 
   \sqrt{\left[\BFKpd{x}{\ }\left(\BFKpd{x}{\BFKav \BFKtu}
         +\BFKpd{y}{\BFKav \BFKtv}\right)\right]^2
     +\left[\BFKpd{y}{\ }\left(\BFKpd{x}{\BFKav \BFKtu}
         +\BFKpd{y}{\BFKav \BFKtv}\right)\right]^2}\end{aligned}

 Whether there is any physical rationale for this correction is unclear
at the moment, but the numerical consequences are good. The divergence
in flows with the grid scale larger or comparable to the Rossby radius
is typically much smaller than the vorticity, so this adjustment only
rarely adjusts the viscosity if :math:`{\sf
  viscC2LeithD}={\sf viscC2Leith}`. However, the rare regions where this
viscosity acts are often the locations for the largest vales of vertical
velocity in the domain. Since the CFL condition on vertical velocity is
often what sets the maximum timestep, this viscosity may substantially
increase the allowable timestep without severely compromising the verity
of the simulation. Tests have shown that in some calculations, a
timestep three times larger was allowed when
:math:`{\sf viscC2LeithD}={\sf viscC2Leith}`.

Courant–Freidrichs–Lewy Constraint on Viscosity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whatever viscosities are used in the model, the choice is constrained by
gridscale and timestep by the Courant–Freidrichs–Lewy (CFL) constraint
on stability:

.. math::

   \begin{aligned}
   A_h & < & \frac{L^2}{4\Delta t} \\
   A_4 & \le & \frac{L^4}{32\Delta t}\end{aligned}

 The viscosities may be automatically limited to be no greater than
these values in MITgcm by specifying viscAhGridMax\ :math:`<1` and
viscA4GridMax\ :math:`<1`. Similarly-scaled minimum values of
viscosities are provided by viscAhGridMin and viscA4GridMin, which if
used, should be set to values :math:`\ll 1`. :math:`L` is roughly the
gridscale (see below).

Following :raw-latex:`\cite{griffies:00}`, we note that there is a
factor of :math:`\Delta
x^2/8` difference between the harmonic and biharmonic viscosities. Thus,
whenever a non-dimensional harmonic coefficient is used in the MITgcm
(*eg.* viscAhGridMax\ :math:`<1`), the biharmonic equivalent is scaled
so that the same non-dimensional value can be used (*eg.*
viscA4GridMax\ :math:`<1`).

Biharmonic Viscosity
~~~~~~~~~~~~~~~~~~~~

:raw-latex:`\cite{ho78}` suggested that eddy viscosities ought to be
focuses on the dynamics at the grid scale, as larger motions would be
’resolved’. To enhance the scale selectivity of the viscous operator, he
suggested a biharmonic eddy viscosity instead of a harmonic (or
Laplacian) viscosity:

.. math::

   \begin{aligned}
   \left({\BFKav{\BFKDt \BFKtu}}-{\BFKaDt \BFKatu}\right)\approx
   \frac{-\nabla^4_h{\BFKatu}}{\BFKRe_4}
   +\frac{\BFKpds{z}{\BFKatu}}{\BFKRe_v}\label{eq:bieddyvisc}, & &
   \left({\BFKav{\BFKDt \BFKtv}}-{\BFKaDt \BFKatv}\right)\approx
   \frac{-\nabla^4_h{\BFKatv}}{\BFKRe_4}
   +\frac{\BFKpds{z}{\BFKatv}}{\BFKRe_v}\nonumber\\
   \left(\BFKav{\BFKDt w}-\BFKaDt
     {\BFKav{w}}\right)\approx\frac{-\nabla^4_h\BFKav
     w}{\BFKRe_4}+\frac{\BFKpds{z}{\BFKav w}}{\BFKRe_v}\nonumber, & &
   \left(\BFKav{\BFKDt{b}}-\BFKaDt{\ \BFKav b} \right)\approx
   \frac{-\nabla^4_h \BFKav b}{\Pr\BFKRe_4}
   +\frac{\BFKpds{z} {\BFKav b}}{\Pr\BFKRe_v}\nonumber\end{aligned}

 :raw-latex:`\cite{griffies:00}` propose that if one scales the
biharmonic viscosity by stability considerations, then the biharmonic
viscous terms will be similarly active to harmonic viscous terms at the
gridscale of the model, but much less active on larger scale motions.
Similarly, a biharmonic diffusivity can be used for less diffusive
flows.

In practice, biharmonic viscosity and diffusivity allow a less viscous,
yet numerically stable, simulation than harmonic viscosity and
diffusivity. However, there is no physical rationale for such operators
being of leading order, and more boundary conditions must be specified
than for the harmonic operators. If one considers the approximations of
[eq:eddyvisc] and [eq:bieddyvisc] to be terms in the Taylor series
expansions of the eddy terms as functions of the large-scale gradient,
then one can argue that both harmonic and biharmonic terms would occur
in the series, and the only question is the choice of coefficients.
Using biharmonic viscosity alone implies that one zeros the first
non-vanishing term in the Taylor series, which is unsupported by any
fluid theory or observation.

Nonetheless, MITgcm supports a plethora of biharmonic viscosities and
diffusivities, which are controlled with parameters named similarly to
the harmonic viscosities and diffusivities with the substitution
:math:`h\rightarrow 4`. MITgcm also supports a biharmonic Leith and
Smagorinsky viscosities:

.. math::

   \begin{aligned}
   A_{4Smag} & = & 
   \left(\frac{{\sf viscC4Smag}}{\pi}\right)^2\frac{L^4}{8}|D| \\
   A_{4Leith} & = & 
   \frac{L^5}{8}\sqrt{\left(\frac{{\sf viscC4Leith}}{\pi}\right)^6
     |\nabla \BFKav \omega_3|^2
     +\left(\frac{{\sf viscC4LeithD}}{\pi}\right)^6
     |\nabla \nabla\cdot \BFKav {\bf \BFKtu}_h|^2}\end{aligned}

 However, it should be noted that unlike the harmonic forms, the
biharmonic scaling does not easily relate to whether energy-dissipation
or enstrophy-dissipation scales are resolved. If similar arguments are
used to estimate these scales and scale them to the gridscale, the
resulting biharmonic viscosities should be:

.. math::

   \begin{aligned}
   A_{4Smag} & = & 
   \left(\frac{{\sf viscC4Smag}}{\pi}\right)^5L^5
   |\nabla^2\BFKav {\bf \BFKtu}_h| \\
   A_{4Leith} & = & 
   L^6\sqrt{\left(\frac{{\sf viscC4Leith}}{\pi}\right)^{12}
     |\nabla^2 \BFKav \omega_3|^2
     +\left(\frac{{\sf viscC4LeithD}}{\pi}\right)^{12}
     |\nabla^2 \nabla\cdot \BFKav {\bf \BFKtu}_h|^2}\end{aligned}

 Thus, the biharmonic scaling suggested by
:raw-latex:`\cite{griffies:00}` implies:

.. math::

   \begin{aligned}
   |D| & \propto &  L|\nabla^2\BFKav {\bf \BFKtu}_h|\\
   |\nabla \BFKav \omega_3| & \propto & L|\nabla^2 \BFKav \omega_3|\end{aligned}

 It is not at all clear that these assumptions ought to hold. Only the
:raw-latex:`\cite{griffies:00}` forms are currently implemented in
MITgcm.

Selection of Length Scale
~~~~~~~~~~~~~~~~~~~~~~~~~

Above, the length scale of the grid has been denoted :math:`L`. However,
in strongly anisotropic grids, :math:`L_x` and :math:`L_y` will be quite
different in some locations. In that case, the CFL condition suggests
that the minimum of :math:`L_x` and :math:`L_y` be used. On the other
hand, other viscosities which involve whether a particular wavelength is
’resolved’ might be better suited to use the maximum of :math:`L_x` and
:math:`L_y`. Currently, MITgcm uses useAreaViscLength to select between
two options. If false, the geometric mean of :math:`L^2_x` and
:math:`L^2_y` is used for all viscosities, which is closer to the
minimum and occurs naturally in the CFL constraint. If useAreaViscLength
is true, then the square root of the area of the grid cell is used.

Mercator, Nondimensional Equations
----------------------------------

The rotating, incompressible, Boussinesq equations of motion
:raw-latex:`\cite{Gill1982}` on a sphere can be written in Mercator
projection about a latitude :math:`\theta_0` and geopotential height
:math:`z=r-r_0`. The nondimensional form of these equations is:

.. math::

   \BFKRo\BFKDt\BFKtu- \frac{\BFKtv
     \sin\theta}{\sin\theta_0}+\BFKMr\BFKpd{x}{\pi}
   +\frac{\lambda\BFKFr^2\BFKMr\cos \theta}{\mu\sin\theta_0} w
   = -\frac{\BFKFr^2\BFKMr \BFKtu w}{r/H}
   +\frac{\BFKRo{\bf \hat x}\cdot\nabla^2{\bf u}}{\BFKRe}

.. math::

   \BFKRo\BFKDt\BFKtv+
   \frac{\BFKtu\sin\theta}{\sin\theta_0}+\BFKMr\BFKpd{y}{\pi}
   = -\frac{\mu\BFKRo\tan\theta(\BFKtu^2+\BFKtv^2)}{r/L} 
   -\frac{\BFKFr^2\BFKMr \BFKtv w}{r/H}
   +\frac{\BFKRo{\bf \hat y}\cdot\nabla^2{\bf u}}{\BFKRe}

.. math::

   \begin{aligned}
   \BFKFr^2\lambda^2\BFKDt w -b+\BFKpd{z}{\pi}
   -\frac{\lambda\cot \theta_0 \BFKtu}{\BFKMr}
   & = & \frac{\lambda\mu^2(\BFKtu^2+\BFKtv^2)}{\BFKMr(r/L)}
   +\frac{\BFKFr^2\lambda^2{\bf \hat z}\cdot\nabla^2{\bf u}}{\BFKRe} \\
   \BFKDt b+w & = & \frac{\nabla^2 b}{\Pr\BFKRe}\nonumber \\
   \mu^2\left(\BFKpd x\BFKtu  + \BFKpd y\BFKtv \right)+\BFKpd z w 
   & = & 0\end{aligned}

 Where

.. math::

   \mu\equiv\frac{\cos\theta_0}{\cos\theta},\ \ \
   \BFKtu=\frac{u^*}{V\mu},\ \ \  \BFKtv=\frac{v^*}{V\mu}

.. math::

   f_0\equiv2\Omega\sin\theta_0,\ \ \  
   %,\ \ \  \BFKDt\  \equiv \mu^2\left(\BFKtu\BFKpd x\  
   %+\BFKtv \BFKpd y\  \right)+\frac{\BFKFr^2\BFKMr}{\BFKRo} w\BFKpd z\  
   \frac{D}{Dt}  \equiv \mu^2\left(\BFKtu\frac{\partial}{\partial x}  
   +\BFKtv \frac{\partial}{\partial y}  \right)
   +\frac{\BFKFr^2\BFKMr}{\BFKRo} w\frac{\partial}{\partial z}

.. math::

   x\equiv \frac{r}{L} \phi \cos \theta_0, \ \ \   
   y\equiv \frac{r}{L} \int_{\theta_0}^\theta
   \frac{\cos \theta_0 \BFKd \theta'}{\cos\theta'}, \ \ \   
   z\equiv \lambda\frac{r-r_0}{L}

.. math:: t^*=t \frac{L}{V},\ \ \  b^*= b\frac{V f_0\BFKMr}{\lambda}

.. math::

   \pi^*=\pi V f_0 L\BFKMr,\ \ \  
   w^*=w V \frac{\BFKFr^2\lambda\BFKMr}{\BFKRo}

.. math:: \BFKRo\equiv\frac{V}{f_0 L},\ \ \  \BFKMr\equiv \max[1,\BFKRo]

.. math::

   \BFKFr\equiv\frac{V}{N \lambda L}, \ \ \   
   \BFKRe\equiv\frac{VL}{\nu}, \ \ \   
   \BFKPr\equiv\frac{\nu}{\kappa}

 Dimensional variables are denoted by an asterisk where necessary. If we
filter over a grid scale typical for ocean models (:math:`1m<L<100km`,
:math:`0.0001<\lambda<1`, :math:`0.001m/s <V<1 m/s`,
:math:`f_0<0.0001 s^{-1}`, :math:`0.01
s^{-1}<N<0.0001 s^{-1}`), these equations are very well approximated by

.. math::

   \begin{aligned}
   \BFKRo{\BFKDt\BFKtu}- \frac{\BFKtv
     \sin\theta}{\sin\theta_0}+\BFKMr\BFKpd{x}{\pi}
   & =& -\frac{\lambda\BFKFr^2\BFKMr\cos \theta}{\mu\sin\theta_0} w
   +\frac{\BFKRo\nabla^2{\BFKtu}}{\BFKRe} \\
   \BFKRo\BFKDt\BFKtv+
   \frac{\BFKtu\sin\theta}{\sin\theta_0}+\BFKMr\BFKpd{y}{\pi}
   & = & \frac{\BFKRo\nabla^2{\BFKtv}}{\BFKRe} \\
   \BFKFr^2\lambda^2\BFKDt w -b+\BFKpd{z}{\pi}
   & = & \frac{\lambda\cot \theta_0 \BFKtu}{\BFKMr}
   +\frac{\BFKFr^2\lambda^2\nabla^2w}{\BFKRe} \\
   \BFKDt b+w & = & \frac{\nabla^2 b}{\Pr\BFKRe} \\
   \mu^2\left(\BFKpd x\BFKtu + \BFKpd y\BFKtv \right)+\BFKpd z w
   & = & 0 \\
   \nabla^2 & \approx & \left(\frac{\partial^2}{\partial x^2}
     +\frac{\partial^2}{\partial y^2}
     +\frac{\partial^2}{\lambda^2\partial z^2}\right)\end{aligned}

 Neglecting the non-frictional terms on the right-hand side is usually
called the ’traditional’ approximation. It is appropriate, with either
large aspect ratio or far from the tropics. This approximation is used
here, as it does not affect the form of the eddy stresses which is the
main topic. The frictional terms are preserved in this approximate form
for later comparison with eddy stresses.
