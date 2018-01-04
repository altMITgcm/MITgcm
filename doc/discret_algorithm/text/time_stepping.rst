This chapter lays out the numerical schemes that are employed in the
core MITgcm algorithm. Whenever possible links are made to actual
program code in the MITgcm implementation. The chapter begins with a
discussion of the temporal discretization used in MITgcm. This
discussion is followed by sections that describe the spatial
discretization. The schemes employed for momentum terms are described
first, afterwards the schemes that apply to passive and dynamically
active tracers are described.

Time-stepping
=============

<!– CMIREDIR:time-stepping: –>

The equations of motion integrated by the model involve four prognostic
equations for flow, :math:`u` and :math:`v`, temperature,
:math:`\theta`, and salt/moisture, :math:`S`, and three diagnostic
equations for vertical flow, :math:`w`, density/buoyancy,
:math:`\rho`/:math:`b`, and pressure/geo-potential, :math:`\phi_{hyd}`.
In addition, the surface pressure or height may by described by either a
prognostic or diagnostic equation and if non-hydrostatics terms are
included then a diagnostic equation for non-hydrostatic pressure is also
solved. The combination of prognostic and diagnostic equations requires
a model algorithm that can march forward prognostic variables while
satisfying constraints imposed by diagnostic equations.

Since the model comes in several flavors and formulation, it would be
confusing to present the model algorithm exactly as written into code
along with all the switches and optional terms. Instead, we present the
algorithm for each of the basic formulations which are:

#. the semi-implicit pressure method for hydrostatic equations with a
   rigid-lid, variables co-located in time and with Adams-Bashforth
   time-stepping, [it:a]

#. as [it:a]. but with an implicit linear free-surface, [it:b]

#. as [it:a]. or [it:b]. but with variables staggered in time, [it:c]

#. as [it:a]. or [it:b]. but with non-hydrostatic terms included,

#. as [it:b]. or [it:c]. but with non-linear free-surface.

In all the above configurations it is also possible to substitute the
Adams-Bashforth with an alternative time-stepping scheme for terms
evaluated explicitly in time. Since the over-arching algorithm is
independent of the particular time-stepping scheme chosen we will
describe first the over-arching algorithm, known as the pressure method,
with a rigid-lid model in section [sec:pressure-method-rigid-lid]. This
algorithm is essentially unchanged, apart for some coefficients, when
the rigid lid assumption is replaced with a linearized implicit
free-surface, described in section
[sec:pressure-method-linear-backward]. These two flavors of the
pressure-method encompass all formulations of the model as it exists
today. The integration of explicit in time terms is out-lined in section
[sec:adams-bashforth] and put into the context of the overall algorithm
in sections [sec:adams-bashforth-sync] and
[sec:adams-bashforth-staggered]. Inclusion of non-hydrostatic terms
requires applying the pressure method in three dimensions instead of two
and this algorithm modification is described in section
[sec:non-hydrostatic]. Finally, the free-surface equation may be treated
more exactly, including non-linear terms, and this is described in
section [sec:nonlinear-freesurface].

Pressure method with rigid-lid
==============================

<!– CMIREDIR:pressure\_method\_rigid\_lid: –>

The horizontal momentum and continuity equations for the ocean
([eq:ocean-mom] and [eq:ocean-cont]), or for the atmosphere
([eq:atmos-mom] and [eq:atmos-cont]), can be summarized by:

.. math::

   \begin{aligned}
   \partial_t u + g \partial_x \eta & = & G_u \\
   \partial_t v + g \partial_y \eta & = & G_v \\
   \partial_x u + \partial_y v + \partial_z w & = & 0\end{aligned}

 where we are adopting the oceanic notation for brevity. All terms in
the momentum equations, except for surface pressure gradient, are
encapsulated in the :math:`G` vector. The continuity equation, when
integrated over the fluid depth, :math:`H`, and with the rigid-lid/no
normal flow boundary conditions applied, becomes:

.. math::

   \partial_x H \widehat{u} + \partial_y H \widehat{v} = 0
   \label{eq:rigid-lid-continuity}

 Here, :math:`H\widehat{u} = \int_H u dz` is the depth integral of
:math:`u`, similarly for :math:`H\widehat{v}`. The rigid-lid
approximation sets :math:`w=0` at the lid so that it does not move but
allows a pressure to be exerted on the fluid by the lid. The horizontal
momentum equations and vertically integrated continuity equation are be
discretized in time and space as follows:

.. math::

   \begin{aligned}
   u^{n+1} + \Delta t g \partial_x \eta^{n+1}
   & = & u^{n} + \Delta t G_u^{(n+1/2)}
   \label{eq:discrete-time-u}
   \\
   v^{n+1} + \Delta t g \partial_y \eta^{n+1}
   & = & v^{n} + \Delta t G_v^{(n+1/2)}
   \label{eq:discrete-time-v}
   \\
     \partial_x H \widehat{u^{n+1}}
   + \partial_y H \widehat{v^{n+1}} & = & 0
   \label{eq:discrete-time-cont-rigid-lid}\end{aligned}

 As written here, terms on the LHS all involve time level :math:`n+1`
and are referred to as implicit; the implicit backward time stepping
scheme is being used. All other terms in the RHS are explicit in time.
The thermodynamic quantities are integrated forward in time in parallel
with the flow and will be discussed later. For the purposes of
describing the pressure method it suffices to say that the hydrostatic
pressure gradient is explicit and so can be included in the vector
:math:`G`.

Substituting the two momentum equations into the depth integrated
continuity equation eliminates :math:`u^{n+1}` and :math:`v^{n+1}`
yielding an elliptic equation for :math:`\eta^{n+1}`. Equations
[eq:discrete-time-u], [eq:discrete-time-v] and
[eq:discrete-time-cont-rigid-lid] can then be re-arranged as follows:

.. math::

   \begin{aligned}
   u^{*} & = & u^{n} + \Delta t G_u^{(n+1/2)} \label{eq:ustar-rigid-lid} \\
   v^{*} & = & v^{n} + \Delta t G_v^{(n+1/2)} \label{eq:vstar-rigid-lid} \\
     \partial_x \Delta t g H \partial_x \eta^{n+1}
   + \partial_y \Delta t g H \partial_y \eta^{n+1}
   & = &
     \partial_x H \widehat{u^{*}}
   + \partial_y H \widehat{v^{*}} \label{eq:elliptic}
   \\
   u^{n+1} & = & u^{*} - \Delta t g \partial_x \eta^{n+1} \label{eq:un+1-rigid-lid}\\
   v^{n+1} & = & v^{*} - \Delta t g \partial_y \eta^{n+1} \label{eq:vn+1-rigid-lid}\end{aligned}

 Equations [eq:ustar-rigid-lid] to [eq:vn+1-rigid-lid], solved
sequentially, represent the pressure method algorithm used in the model.
The essence of the pressure method lies in the fact that any explicit
prediction for the flow would lead to a divergence flow field so a
pressure field must be found that keeps the flow non-divergent over each
step of the integration. The particular location in time of the pressure
field is somewhat ambiguous; in Fig. [fig:pressure-method-rigid-lid] we
depicted as co-located with the future flow field (time level
:math:`n+1`) but it could equally have been drawn as staggered in time
with the flow.

The correspondence to the code is as follows:

-  the prognostic phase, equations [eq:ustar-rigid-lid] and
   [eq:vstar-rigid-lid], stepping forward :math:`u^n` and :math:`v^n` to
   :math:`u^{*}` and :math:`v^{*}` is coded in

-  the vertical integration, :math:`H \widehat{u^*}` and :math:`H
   \widehat{v^*}`, divergence and inversion of the elliptic operator in
   equation [eq:elliptic] is coded in

-  finally, the new flow field at time level :math:`n+1` given by
   equations [eq:un+1-rigid-lid] and [eq:vn+1-rigid-lid] is calculated
   in .

The calling tree for these routines is given in
Fig. [fig:call-tree-pressure-method].

In general, the horizontal momentum time-stepping can contain some terms
that are treated implicitly in time, such as the vertical viscosity when
using the backward time-stepping scheme (). The method used to solve
those implicit terms is provided in section
[sec:implicit-backward-stepping], and modifies equations
[eq:discrete-time-u] and [eq:discrete-time-v] to give:

.. math::

   \begin{aligned}
   u^{n+1} - \Delta t \partial_z A_v \partial_z u^{n+1}
   + \Delta t g \partial_x \eta^{n+1} & = & u^{n} + \Delta t G_u^{(n+1/2)}
   \\
   v^{n+1} - \Delta t \partial_z A_v \partial_z v^{n+1}
   + \Delta t g \partial_y \eta^{n+1} & = & v^{n} + \Delta t G_v^{(n+1/2)}\end{aligned}

Pressure method with implicit linear free-surface
=================================================

<!– CMIREDIR:pressure\_method\_linear\_backward: –>

The rigid-lid approximation filters out external gravity waves
subsequently modifying the dispersion relation of barotropic Rossby
waves. The discrete form of the elliptic equation has some zero
eigen-values which makes it a potentially tricky or inefficient problem
to solve.

The rigid-lid approximation can be easily replaced by a linearization of
the free-surface equation which can be written:

.. math::

   \partial_t \eta + \partial_x H \widehat{u} + \partial_y H \widehat{v} = P-E+R
   \label{eq:linear-free-surface=P-E}

 which differs from the depth integrated continuity equation with
rigid-lid ([eq:rigid-lid-continuity]) by the time-dependent term and
fresh-water source term.

Equation [eq:discrete-time-cont-rigid-lid] in the rigid-lid pressure
method is then replaced by the time discretization of
[eq:linear-free-surface=P-E] which is:

.. math::

   \eta^{n+1}
   + \Delta t \partial_x H \widehat{u^{n+1}}
   + \Delta t \partial_y H \widehat{v^{n+1}}
   =
   \eta^{n}
   + \Delta t ( P - E )
   \label{eq:discrete-time-backward-free-surface}

 where the use of flow at time level :math:`n+1` makes the method
implicit and backward in time. This is the preferred scheme since it
still filters the fast unresolved wave motions by damping them. A
centered scheme, such as Crank-Nicholson (see section
[sec:freesurf-CrankNick]), would alias the energy of the fast modes onto
slower modes of motion.

As for the rigid-lid pressure method, equations [eq:discrete-time-u],
[eq:discrete-time-v] and [eq:discrete-time-backward-free-surface] can be
re-arranged as follows:

.. math::

   \begin{aligned}
   u^{*} & = & u^{n} + \Delta t G_u^{(n+1/2)} \label{eq:ustar-backward-free-surface} \\
   v^{*} & = & v^{n} + \Delta t G_v^{(n+1/2)} \label{eq:vstar-backward-free-surface} \\
   \eta^* & = & \epsilon_{fs} \left( \eta^{n} + \Delta t (P-E) \right)
            - \Delta t \left( \partial_x H \widehat{u^{*}}
                            + \partial_y H \widehat{v^{*}} \right)
   \\
     \partial_x g H \partial_x \eta^{n+1}
   & + & \partial_y g H \partial_y \eta^{n+1}
    - \frac{\epsilon_{fs} \eta^{n+1}}{\Delta t^2}
    =
   - \frac{\eta^*}{\Delta t^2}
   \label{eq:elliptic-backward-free-surface}
   \\
   u^{n+1} & = & u^{*} - \Delta t g \partial_x \eta^{n+1} \label{eq:un+1-backward-free-surface}\\
   v^{n+1} & = & v^{*} - \Delta t g \partial_y \eta^{n+1} \label{eq:vn+1-backward-free-surface}\end{aligned}

 Equations [eq:ustar-backward-free-surface]
to [eq:vn+1-backward-free-surface], solved sequentially, represent the
pressure method algorithm with a backward implicit, linearized free
surface. The method is still formerly a pressure method because in the
limit of large :math:`\Delta t` the rigid-lid method is recovered.
However, the implicit treatment of the free-surface allows the flow to
be divergent and for the surface pressure/elevation to respond on a
finite time-scale (as opposed to instantly). To recover the rigid-lid
formulation, we introduced a switch-like parameter,
:math:`\epsilon_{fs}` (), which selects between the free-surface and
rigid-lid; :math:`\epsilon_{fs}=1` allows the free-surface to evolve;
:math:`\epsilon_{fs}=0` imposes the rigid-lid. The evolution in time and
location of variables is exactly as it was for the rigid-lid model so
that Fig. [fig:pressure-method-rigid-lid] is still applicable.
Similarly, the calling sequence, given in
Fig. [fig:call-tree-pressure-method], is as for the pressure-method.

Explicit time-stepping: Adams-Bashforth
=======================================

<!– CMIREDIR:adams\_bashforth: –>

In describing the the pressure method above we deferred describing the
time discretization of the explicit terms. We have historically used the
quasi-second order Adams-Bashforth method for all explicit terms in both
the momentum and tracer equations. This is still the default mode of
operation but it is now possible to use alternate schemes for tracers
(see section [sec:tracer-advection]).

In the previous sections, we summarized an explicit scheme as:

.. math::

   \tau^{*} = \tau^{n} + \Delta t G_\tau^{(n+1/2)}
   \label{eq:taustar}

 where :math:`\tau` could be any prognostic variable (:math:`u`,
:math:`v`, :math:`\theta` or :math:`S`) and :math:`\tau^*` is an
explicit estimate of :math:`\tau^{n+1}` and would be exact if not for
implicit-in-time terms. The parenthesis about :math:`n+1/2` indicates
that the term is explicit and extrapolated forward in time and for this
we use the quasi-second order Adams-Bashforth method:

.. math::

   G_\tau^{(n+1/2)} = ( 3/2 + \epsilon_{AB}) G_\tau^n
   - ( 1/2 + \epsilon_{AB}) G_\tau^{n-1}
   \label{eq:adams-bashforth2}

 This is a linear extrapolation, forward in time, to
:math:`t=(n+1/2+{\epsilon_{AB}})\Delta t`. An extrapolation to the
mid-point in time, :math:`t=(n+1/2)\Delta t`, corresponding to
:math:`\epsilon_{AB}=0`, would be second order accurate but is weakly
unstable for oscillatory terms. A small but finite value for
:math:`\epsilon_{AB}` stabilizes the method. Strictly speaking, damping
terms such as diffusion and dissipation, and fixed terms (forcing), do
not need to be inside the Adams-Bashforth extrapolation. However, in the
current code, it is simpler to include these terms and this can be
justified if the flow and forcing evolves smoothly. Problems can, and
do, arise when forcing or motions are high frequency and this
corresponds to a reduced stability compared to a simple forward
time-stepping of such terms. The model offers the possibility to leave
the tracer and momentum forcing terms and the dissipation terms outside
the Adams-Bashforth extrapolation, by turning off the logical flags
(parameter file *data*, namelist *PARM01*, default value = True). (,
default=0, , default=0) and (parameter file *data*, namelist *PARM01*,
default value = True). respectively.

A stability analysis for an oscillation equation should be given at this
point.

A stability analysis for a relaxation equation should be given at this
point.

Implicit time-stepping: backward method
=======================================

<!– CMIREDIR:implicit\_time-stepping\_backward: –>

Vertical diffusion and viscosity can be treated implicitly in time using
the backward method which is an intrinsic scheme. Recently, the option
to treat the vertical advection implicitly has been added, but not yet
tested; therefore, the description hereafter is limited to diffusion and
viscosity. For tracers, the time discretized equation is:

.. math::

   \tau^{n+1} - \Delta t \partial_r \kappa_v \partial_r \tau^{n+1} =
   \tau^{n} + \Delta t G_\tau^{(n+1/2)}
   \label{eq:implicit-diffusion}

 where :math:`G_\tau^{(n+1/2)}` is the remaining explicit terms
extrapolated using the Adams-Bashforth method as described above.
Equation [eq:implicit-diffusion] can be split split into:

.. math::

   \begin{aligned}
   \tau^* & = & \tau^{n} + \Delta t G_\tau^{(n+1/2)}
   \label{eq:taustar-implicit} \\
   \tau^{n+1} & = & {\cal L}_\tau^{-1} ( \tau^* )
   \label{eq:tau-n+1-implicit}\end{aligned}

 where :math:`{\cal L}_\tau^{-1}` is the inverse of the operator

.. math:: {\cal L}_\tau = \left[ 1 + \Delta t \partial_r \kappa_v \partial_r \right]

 Equation [eq:taustar-implicit] looks exactly as [eq:taustar] while
[eq:tau-n+1-implicit] involves an operator or matrix inversion. By
re-arranging [eq:implicit-diffusion] in this way we have cast the method
as an explicit prediction step and an implicit step allowing the latter
to be inserted into the over all algorithm with minimal interference.

Fig. [fig:call-tree-adams-bashforth] shows the calling sequence for
stepping forward a tracer variable such as temperature.

In order to fit within the pressure method, the implicit viscosity must
not alter the barotropic flow. In other words, it can only redistribute
momentum in the vertical. The upshot of this is that although vertical
viscosity may be backward implicit and unconditionally stable, no-slip
boundary conditions may not be made implicit and are thus cast as a an
explicit drag term.

Synchronous time-stepping: variables co-located in time
=======================================================

<!– CMIREDIR:adams\_bashforth\_sync: –>

The Adams-Bashforth extrapolation of explicit tendencies fits neatly
into the pressure method algorithm when all state variables are
co-located in time. Fig. [fig:adams-bashforth-sync] illustrates the
location of variables in time and the evolution of the algorithm with
time. The algorithm can be represented by the sequential solution of the
follow equations:

.. math::

   \begin{aligned}
   G_{\theta,S}^{n} & = & G_{\theta,S} ( u^n, \theta^n, S^n )
   \label{eq:Gt-n-sync} \\
   G_{\theta,S}^{(n+1/2)} & = & (3/2+\epsilon_{AB}) G_{\theta,S}^{n}-(1/2+\epsilon_{AB}) G_{\theta,S}^{n-1}
   \label{eq:Gt-n+5-sync} \\
   (\theta^*,S^*) & = & (\theta^{n},S^{n}) + \Delta t G_{\theta,S}^{(n+1/2)}
   \label{eq:tstar-sync} \\
   (\theta^{n+1},S^{n+1}) & = & {\cal L}^{-1}_{\theta,S} (\theta^*,S^*)
   \label{eq:t-n+1-sync} \\
   \phi^n_{hyd} & = & \int b(\theta^n,S^n) dr
   \label{eq:phi-hyd-sync} \\
   \vec{\bf G}_{\vec{\bf v}}^{n} & = & \vec{\bf G}_{\vec{\bf v}} ( \vec{\bf v}^n, \phi^n_{hyd} )
   \label{eq:Gv-n-sync} \\
   \vec{\bf G}_{\vec{\bf v}}^{(n+1/2)} & = & (3/2 + \epsilon_{AB} ) \vec{\bf G}_{\vec{\bf v}}^{n} - (1/2 + \epsilon_{AB} ) \vec{\bf G}_{\vec{\bf v}}^{n-1}
   \label{eq:Gv-n+5-sync} \\
   \vec{\bf v}^{*} & = & \vec{\bf v}^{n} + \Delta t \vec{\bf G}_{\vec{\bf v}}^{(n+1/2)}
   \label{eq:vstar-sync} \\
   \vec{\bf v}^{**} & = & {\cal L}_{\vec{\bf v}}^{-1} ( \vec{\bf v}^* )
   \label{eq:vstarstar-sync} \\
   \eta^* & = & \epsilon_{fs} \left( \eta^{n} + \Delta t (P-E) \right)- \Delta t
     \nabla \cdot H \widehat{ \vec{\bf v}^{**} }
   \label{eq:nstar-sync} \\
   \nabla \cdot g H \nabla \eta^{n+1} & - & \frac{\epsilon_{fs} \eta^{n+1}}{\Delta t^2}
   ~ = ~ - \frac{\eta^*}{\Delta t^2}
   \label{eq:elliptic-sync} \\
   \vec{\bf v}^{n+1} & = & \vec{\bf v}^{**} - \Delta t g \nabla \eta^{n+1}
   \label{eq:v-n+1-sync}\end{aligned}

 Fig. [fig:adams-bashforth-sync] illustrates the location of variables
in time and evolution of the algorithm with time. The Adams-Bashforth
extrapolation of the tracer tendencies is illustrated by the dashed
arrow, the prediction at :math:`n+1` is indicated by the solid arc.
Inversion of the implicit terms, :math:`{\cal
L}^{-1}_{\theta,S}`, then yields the new tracer fields at :math:`n+1`.
All these operations are carried out in subroutine *THERMODYNAMICS* an
subsidiaries, which correspond to equations [eq:Gt-n-sync] to
[eq:t-n+1-sync]. Similarly illustrated is the Adams-Bashforth
extrapolation of accelerations, stepping forward and solving of implicit
viscosity and surface pressure gradient terms, corresponding to
equations [eq:Gv-n-sync] to [eq:v-n+1-sync]. These operations are
carried out in subroutines *DYNAMCIS*, *SOLVE\_FOR\_PRESSURE* and
*MOMENTUM\_CORRECTION\_STEP*. This, then, represents an entire algorithm
for stepping forward the model one time-step. The corresponding calling
tree is given in [fig:call-tree-adams-bashforth-sync].

Staggered baroclinic time-stepping
==================================

<!– CMIREDIR:adams\_bashforth\_staggered: –>

For well stratified problems, internal gravity waves may be the limiting
process for determining a stable time-step. In the circumstance, it is
more efficient to stagger in time the thermodynamic variables with the
flow variables. Fig. [fig:adams-bashforth-staggered] illustrates the
staggering and algorithm. The key difference between this and
Fig. [fig:adams-bashforth-sync] is that the thermodynamic variables are
solved after the dynamics, using the recently updated flow field. This
essentially allows the gravity wave terms to leap-frog in time giving
second order accuracy and more stability.

The essential change in the staggered algorithm is that the
thermodynamics solver is delayed from half a time step, allowing the use
of the most recent velocities to compute the advection terms. Once the
thermodynamics fields are updated, the hydrostatic pressure is computed
to step forward the dynamics. Note that the pressure gradient must also
be taken out of the Adams-Bashforth extrapolation. Also, retaining the
integer time-levels, :math:`n` and :math:`n+1`, does not give a user the
sense of where variables are located in time. Instead, we re-write the
entire algorithm, [eq:Gt-n-sync] to [eq:v-n+1-sync], annotating the
position in time of variables appropriately:

.. math::

   \begin{aligned}
   \phi^{n}_{hyd} & = & \int b(\theta^{n},S^{n}) dr
   \label{eq:phi-hyd-staggered} \\
   \vec{\bf G}_{\vec{\bf v}}^{n-1/2} & = & \vec{\bf G}_{\vec{\bf v}} ( \vec{\bf v}^{n-1/2} )
   \label{eq:Gv-n-staggered} \\
   \vec{\bf G}_{\vec{\bf v}}^{(n)} & = & (3/2 + \epsilon_{AB} ) \vec{\bf G}_{\vec{\bf v}}^{n-1/2} - (1/2 + \epsilon_{AB} ) \vec{\bf G}_{\vec{\bf v}}^{n-3/2}
   \label{eq:Gv-n+5-staggered} \\
   \vec{\bf v}^{*} & = & \vec{\bf v}^{n-1/2} + \Delta t \left( \vec{\bf G}_{\vec{\bf v}}^{(n)} - \nabla \phi_{hyd}^{n} \right)
   \label{eq:vstar-staggered} \\
   \vec{\bf v}^{**} & = & {\cal L}_{\vec{\bf v}}^{-1} ( \vec{\bf v}^* )
   \label{eq:vstarstar-staggered} \\
   \eta^* & = & \epsilon_{fs} \left( \eta^{n-1/2} + \Delta t (P-E)^n \right)- \Delta t
     \nabla \cdot H \widehat{ \vec{\bf v}^{**} }
   \label{eq:nstar-staggered} \\
   \nabla \cdot g H \nabla \eta^{n+1/2} & - & \frac{\epsilon_{fs} \eta^{n+1/2}}{\Delta t^2}
   ~ = ~ - \frac{\eta^*}{\Delta t^2}
   \label{eq:elliptic-staggered} \\
   \vec{\bf v}^{n+1/2} & = & \vec{\bf v}^{**} - \Delta t g \nabla \eta^{n+1/2}
   \label{eq:v-n+1-staggered} \\
   G_{\theta,S}^{n} & = & G_{\theta,S} ( u^{n+1/2}, \theta^{n}, S^{n} )
   \label{eq:Gt-n-staggered} \\
   G_{\theta,S}^{(n+1/2)} & = & (3/2+\epsilon_{AB}) G_{\theta,S}^{n}-(1/2+\epsilon_{AB}) G_{\theta,S}^{n-1}
   \label{eq:Gt-n+5-staggered} \\
   (\theta^*,S^*) & = & (\theta^{n},S^{n}) + \Delta t G_{\theta,S}^{(n+1/2)}
   \label{eq:tstar-staggered} \\
   (\theta^{n+1},S^{n+1}) & = & {\cal L}^{-1}_{\theta,S} (\theta^*,S^*)
   \label{eq:t-n+1-staggered}\end{aligned}

 The corresponding calling tree is given in
[fig:call-tree-adams-bashforth-staggered]. The staggered algorithm is
activated with the run-time flag **staggerTimeStep**\ *=.TRUE.* in
parameter file *data*, namelist *PARM01*.

The only difficulty with this approach is apparent in equation
[eq:Gt-n-staggered] and illustrated by the dotted arrow connecting
:math:`u,v^{n+1/2}` with :math:`G_\theta^{n}`. The flow used to advect
tracers around is not naturally located in time. This could be avoided
by applying the Adams-Bashforth extrapolation to the tracer field itself
and advecting that around but this approach is not yet available. We’re
not aware of any detrimental effect of this feature. The difficulty lies
mainly in interpretation of what time-level variables and terms
correspond to.

Non-hydrostatic formulation
===========================

<!– CMIREDIR:non-hydrostatic\_formulation: –>

The non-hydrostatic formulation re-introduces the full vertical momentum
equation and requires the solution of a 3-D elliptic equations for
non-hydrostatic pressure perturbation. We still integrate vertically for
the hydrostatic pressure and solve a 2-D elliptic equation for the
surface pressure/elevation for this reduces the amount of work needed to
solve for the non-hydrostatic pressure.

The momentum equations are discretized in time as follows:

.. math::

   \begin{aligned}
   \frac{1}{\Delta t} u^{n+1} + g \partial_x \eta^{n+1} + \partial_x \phi_{nh}^{n+1}
   & = & \frac{1}{\Delta t} u^{n} + G_u^{(n+1/2)} \label{eq:discrete-time-u-nh} \\
   \frac{1}{\Delta t} v^{n+1} + g \partial_y \eta^{n+1} + \partial_y \phi_{nh}^{n+1}
   & = & \frac{1}{\Delta t} v^{n} + G_v^{(n+1/2)} \label{eq:discrete-time-v-nh} \\
   \frac{1}{\Delta t} w^{n+1} + \partial_r \phi_{nh}^{n+1}
   & = & \frac{1}{\Delta t} w^{n} + G_w^{(n+1/2)} \label{eq:discrete-time-w-nh}\end{aligned}

 which must satisfy the discrete-in-time depth integrated continuity,
equation [eq:discrete-time-backward-free-surface] and the local
continuity equation

.. math::

   \partial_x u^{n+1} + \partial_y v^{n+1} + \partial_r w^{n+1} = 0
   \label{eq:non-divergence-nh}

 As before, the explicit predictions for momentum are consolidated as:

.. math::

   \begin{aligned}
   u^* & = & u^n + \Delta t G_u^{(n+1/2)} \\
   v^* & = & v^n + \Delta t G_v^{(n+1/2)} \\
   w^* & = & w^n + \Delta t G_w^{(n+1/2)}\end{aligned}

 but this time we introduce an intermediate step by splitting the
tendancy of the flow as follows:

.. math::

   \begin{aligned}
   u^{n+1} = u^{**} - \Delta t \partial_x \phi_{nh}^{n+1}
   & &
   u^{**} = u^{*} - \Delta t g \partial_x \eta^{n+1} \\
   v^{n+1} = v^{**} - \Delta t \partial_y \phi_{nh}^{n+1}
   & &
   v^{**} = v^{*} - \Delta t g \partial_y \eta^{n+1}\end{aligned}

 Substituting into the depth integrated continuity
(equation [eq:discrete-time-backward-free-surface]) gives

.. math::

   \partial_x H \partial_x \left( g \eta^{n+1} + \widehat{\phi}_{nh}^{n+1} \right)
   +
   \partial_y H \partial_y \left( g \eta^{n+1} + \widehat{\phi}_{nh}^{n+1} \right)
    - \frac{\epsilon_{fs}\eta^{n+1}}{\Delta t^2}
   = - \frac{\eta^*}{\Delta t^2}

 which is approximated by equation [eq:elliptic-backward-free-surface]
on the basis that i) :math:`\phi_{nh}^{n+1}` is not yet known and ii)
:math:`\nabla \widehat{\phi}_{nh}
<< g \nabla \eta`. If [eq:elliptic-backward-free-surface] is solved
accurately then the implication is that :math:`\widehat{\phi}_{nh}
\approx 0` so that the non-hydrostatic pressure field does not drive
barotropic motion.

The flow must satisfy non-divergence (equation [eq:non-divergence-nh])
locally, as well as depth integrated, and this constraint is used to
form a 3-D elliptic equations for :math:`\phi_{nh}^{n+1}`:

.. math::

   \partial_{xx} \phi_{nh}^{n+1} + \partial_{yy} \phi_{nh}^{n+1} +
   \partial_{rr} \phi_{nh}^{n+1} =
   \partial_x u^{**} + \partial_y v^{**} + \partial_r w^{*}

The entire algorithm can be summarized as the sequential solution of the
following equations:

.. math::

   \begin{aligned}
   u^{*} & = & u^{n} + \Delta t G_u^{(n+1/2)} \label{eq:ustar-nh} \\
   v^{*} & = & v^{n} + \Delta t G_v^{(n+1/2)} \label{eq:vstar-nh} \\
   w^{*} & = & w^{n} + \Delta t G_w^{(n+1/2)} \label{eq:wstar-nh} \\
   \eta^* ~ = ~ \epsilon_{fs} \left( \eta^{n} + \Delta t (P-E) \right)
   & - & \Delta t \left( \partial_x H \widehat{u^{*}}
                       + \partial_y H \widehat{v^{*}} \right)
   \\
     \partial_x g H \partial_x \eta^{n+1}
   + \partial_y g H \partial_y \eta^{n+1}
   & - & \frac{\epsilon_{fs} \eta^{n+1}}{\Delta t^2}
   ~ = ~
   - \frac{\eta^*}{\Delta t^2}
   \label{eq:elliptic-nh}
   \\
   u^{**} & = & u^{*} - \Delta t g \partial_x \eta^{n+1} \label{eq:unx-nh}\\
   v^{**} & = & v^{*} - \Delta t g \partial_y \eta^{n+1} \label{eq:vnx-nh}\\
   \partial_{xx} \phi_{nh}^{n+1} + \partial_{yy} \phi_{nh}^{n+1} +
   \partial_{rr} \phi_{nh}^{n+1} & = &
   \partial_x u^{**} + \partial_y v^{**} + \partial_r w^{*}  \label{eq:phi-nh}\\
   u^{n+1} & = & u^{**} - \Delta t \partial_x \phi_{nh}^{n+1} \label{eq:un+1-nh}\\
   v^{n+1} & = & v^{**} - \Delta t \partial_y \phi_{nh}^{n+1} \label{eq:vn+1-nh}\\
   \partial_r w^{n+1} & = & - \partial_x u^{n+1} - \partial_y v^{n+1}\end{aligned}

 where the last equation is solved by vertically integrating for
:math:`w^{n+1}`.

Variants on the Free Surface
============================

We now describe the various formulations of the free-surface that
include non-linear forms, implicit in time using Crank-Nicholson,
explicit and [one day] split-explicit. First, we’ll reiterate the
underlying algorithm but this time using the notation consistent with
the more general vertical coordinate :math:`r`. The elliptic equation
for free-surface coordinate (units of :math:`r`), corresponding to
[eq:discrete-time-backward-free-surface], and assuming no
non-hydrostatic effects (:math:`\epsilon_{nh} = 0`) is:

.. math::

   \begin{aligned}
   \epsilon_{fs} {\eta}^{n+1} -
   {\bf \nabla}_h \cdot \Delta t^2 (R_o-R_{fixed}) {\bf \nabla}_h b_s
   {\eta}^{n+1} = {\eta}^*
   \label{eq-solve2D}\end{aligned}

 where

.. math::

   \begin{aligned}
   {\eta}^* = \epsilon_{fs} \: {\eta}^{n} -
   \Delta t {\bf \nabla}_h \cdot \int_{R_{fixed}}^{R_o} \vec{\bf v}^* dr
   \: + \: \epsilon_{fw} \Delta t (P-E)^{n}
   \label{eq-solve2D_rhs}\end{aligned}

Once :math:`{\eta}^{n+1}` has been found, substituting into
[eq:discrete-time-u], [eq:discrete-time-v] yields
:math:`\vec{\bf v}^{n+1}` if the model is hydrostatic
(:math:`\epsilon_{nh}=0`):

.. math::

   \vec{\bf v}^{n+1} = \vec{\bf v}^{*}
   - \Delta t {\bf \nabla}_h b_s {\eta}^{n+1}

This is known as the correction step. However, when the model is
non-hydrostatic (:math:`\epsilon_{nh}=1`) we need an additional step and
an additional equation for :math:`\phi'_{nh}`. This is obtained by
substituting [eq:discrete-time-u-nh], [eq:discrete-time-v-nh] and
[eq:discrete-time-w-nh] into continuity:

.. math::

   \left[ {\bf \nabla}_h^2 + \partial_{rr} \right] {\phi'_{nh}}^{n+1}
   = \frac{1}{\Delta t} \left(
   {\bf \nabla}_h \cdot \vec{\bf v}^{**} + \partial_r \dot{r}^* \right)

 where

.. math:: \vec{\bf v}^{**} = \vec{\bf v}^* - \Delta t {\bf \nabla}_h b_s {\eta}^{n+1}

 Note that :math:`\eta^{n+1}` is also used to update the second RHS term
:math:`\partial_r \dot{r}^* ` since the vertical velocity at the surface
(:math:`\dot{r}_{surf}`) is evaluated as
:math:`(\eta^{n+1} - \eta^n) / \Delta t`.

Finally, the horizontal velocities at the new time level are found by:

.. math::

   \vec{\bf v}^{n+1} = \vec{\bf v}^{**}
   - \epsilon_{nh} \Delta t {\bf \nabla}_h {\phi'_{nh}}^{n+1}

 and the vertical velocity is found by integrating the continuity
equation vertically. Note that, for the convenience of the restart
procedure, the vertical integration of the continuity equation has been
moved to the beginning of the time step (instead of at the end), without
any consequence on the solution.

Regarding the implementation of the surface pressure solver, all
computation are done within the routine *SOLVE\_FOR\_PRESSURE* and its
dependent calls. The standard method to solve the 2D elliptic problem
([eq-solve2D]) uses the conjugate gradient method (routine *CG2D*); the
solver matrix and conjugate gradient operator are only function of the
discretized domain and are therefore evaluated separately, before the
time iteration loop, within *INI\_CG2D*. The computation of the RHS
:math:`\eta^*` is partly done in *CALC\_DIV\_GHAT* and in
*SOLVE\_FOR\_PRESSURE*.

The same method is applied for the non hydrostatic part, using a
conjugate gradient 3D solver (*CG3D*) that is initialized in
*INI\_CG3D*. The RHS terms of 2D and 3D problems are computed together
at the same point in the code.

Crank-Nicolson barotropic time stepping
---------------------------------------

| The full implicit time stepping described previously is
  unconditionally stable but damps the fast gravity waves, resulting in
  a loss of potential energy. The modification presented now allows one
  to combine an implicit part (:math:`\beta,\gamma`) and an explicit
  part (:math:`1-\beta,1-\gamma`) for the surface pressure gradient
  (:math:`\beta`) and for the barotropic flow divergence
  (:math:`\gamma`).
| For instance, :math:`\beta=\gamma=1` is the previous fully implicit
  scheme; :math:`\beta=\gamma=1/2` is the non damping (energy
  conserving), unconditionally stable, Crank-Nicolson scheme;
  :math:`(\beta,\gamma)=(1,0)` or :math:`=(0,1)` corresponds to the
  forward - backward scheme that conserves energy but is only stable for
  small time steps.
| In the code, :math:`\beta,\gamma` are defined as parameters,
  respectively **implicSurfPress**, **implicDiv2DFlow**. They are read
  from the main parameter file “*data*” (namelist *PARM01*) and are set
  by default to 1,1.

| Equations [eq:ustar-backward-free-surface] –
  [eq:vn+1-backward-free-surface] are modified as follows:

  .. math::

     \begin{aligned}
     \frac{ \vec{\bf v}^{n+1} }{ \Delta t }
     + {\bf \nabla}_h b_s [ \beta {\eta}^{n+1} + (1-\beta) {\eta}^{n} ]
     + \epsilon_{nh} {\bf \nabla}_h {\phi'_{nh}}^{n+1}
      = \frac{ \vec{\bf v}^{n} }{ \Delta t }
      + \vec{\bf G}_{\vec{\bf v}} ^{(n+1/2)}
      + {\bf \nabla}_h {\phi'_{hyd}}^{(n+1/2)}\end{aligned}

  .. math::

     \begin{aligned}
     \epsilon_{fs} \frac{ {\eta}^{n+1} - {\eta}^{n} }{ \Delta t}
     + {\bf \nabla}_h \cdot \int_{R_{fixed}}^{R_o}
     [ \gamma \vec{\bf v}^{n+1} + (1-\gamma) \vec{\bf v}^{n}] dr
     = \epsilon_{fw} (P-E)
     \label{eq:eta-n+1-CrankNick}\end{aligned}

   We set

  .. math::

     \begin{aligned}
     \vec{\bf v}^* & = &
     \vec{\bf v} ^{n} + \Delta t \vec{\bf G}_{\vec{\bf v}} ^{(n+1/2)}
     + (\beta-1) \Delta t {\bf \nabla}_h b_s {\eta}^{n}
     + \Delta t {\bf \nabla}_h {\phi'_{hyd}}^{(n+1/2)}
     \\
     {\eta}^* & = &
     \epsilon_{fs} {\eta}^{n} + \epsilon_{fw} \Delta t (P-E)
     - \Delta t {\bf \nabla}_h \cdot \int_{R_{fixed}}^{R_o}
     [ \gamma \vec{\bf v}^* + (1-\gamma) \vec{\bf v}^{n}] dr\end{aligned}
| In the hydrostatic case (:math:`\epsilon_{nh}=0`), allowing us to find
  :math:`{\eta}^{n+1}`, thus:

  .. math::

     \epsilon_{fs} {\eta}^{n+1} -
     {\bf \nabla}_h \cdot \beta\gamma \Delta t^2 b_s (R_o - R_{fixed})
     {\bf \nabla}_h {\eta}^{n+1}
     = {\eta}^*

   and then to compute (*CORRECTION\_STEP*):

  .. math::

     \vec{\bf v}^{n+1} = \vec{\bf v}^{*}
     - \beta \Delta t {\bf \nabla}_h b_s {\eta}^{n+1}

Notes:

#. The RHS term of equation [eq:eta-n+1-CrankNick] corresponds the
   contribution of fresh water flux (P-E) to the free-surface variations
   (:math:`\epsilon_{fw}=1`, **useRealFreshWater**\ *=TRUE* in parameter
   file *data*). In order to remain consistent with the tracer equation,
   specially in the non-linear free-surface formulation, this term is
   also affected by the Crank-Nicolson time stepping. The RHS reads:
   :math:`\epsilon_{fw} ( \gamma (P-E)^{n+1/2} + (1-\gamma) (P-E)^{n-1/2} )`

#. The stability criteria with Crank-Nicolson time stepping for the pure
   linear gravity wave problem in cartesian coordinates is:

   -  :math:`\beta + \gamma < 1` : unstable

   -  :math:`\beta \geq 1/2` and :math:` \gamma \geq 1/2` : stable

   -  :math:`\beta + \gamma \geq 1` : stable if

      .. math:: c_{max}^2 (\beta - 1/2)(\gamma - 1/2) + 1 \geq 0

      .. math::

         \mbox{with }~
         %c^2 = 2 g H {\Delta t}^2
         %(\frac{1-cos 2 \pi / k}{\Delta x^2}
         %+\frac{1-cos 2 \pi / l}{\Delta y^2})
         %

       c\_max = 2 t :math:``

#. A similar mixed forward/backward time-stepping is also available for
   the non-hydrostatic algorithm, with a fraction :math:`\beta_{nh}`
   (:math:` 0 < \beta_{nh} \leq 1`) of the non-hydrostatic pressure
   gradient being evaluated at time step :math:`n+1` (backward in time)
   and the remaining part (:math:`1 - \beta_{nh}`) being evaluated at
   time step :math:`n` (forward in time). The run-time parameter
   **implicitNHPress** corresponding to the implicit fraction
   :math:`\beta_{nh}` of the non-hydrostatic pressure is set by default
   to the implicit fraction :math:`\beta` of surface pressure
   (**implicSurfPress**), but can also be specified independently (in
   main parameter file *data*, namelist *PARM01*).
