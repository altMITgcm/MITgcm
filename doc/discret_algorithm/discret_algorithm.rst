.. _discret_algorithm:

Discretization and Algorithm
****************************

This chapter lays out the numerical schemes that are employed in the
core MITgcm algorithm. Whenever possible links are made to actual
program code in the MITgcm implementation. The chapter begins with a
discussion of the temporal discretization used in MITgcm. This
discussion is followed by sections that describe the spatial
discretization. The schemes employed for momentum terms are described
first, afterwards the schemes that apply to passive and dynamically
active tracers are described.

Notation
========

Because of the particularity of the vertical direction in stratified
fluid context, in this chapter, the vector notations are mostly used for
the horizontal component: the horizontal part of a vector is simply
written :math:`\vec{\bf v}` (instead of :math:`{\bf v_h}` or
:math:`\vec{\mathbf{v}}_{h}` in chapter 1) and a 3D vector is simply
written :math:`\vec{v}` (instead of :math:`\vec{\mathbf{v}}` in chapter
1).

The notations we use to describe the discrete formulation of the model
are summarized as follows.

General notation:

| :math:`\Delta x, \Delta y, \Delta r` grid spacing in X, Y, R directions
|
| :math:`A_c,A_w,A_s,A_{\zeta}` : horizontal area of a grid cell
  surrounding :math:`\theta,u,v,\zeta` point
|
| :math:`{\cal V}_u , {\cal V}_v , {\cal V}_w , {\cal V}_\theta` :
  Volume of the grid box surrounding :math:`u,v,w,\theta` point
|
| :math:`i,j,k` : current index relative to X, Y, R directions
| 
Basic operators:

| :math:`\delta_i` :
  :math:`\delta_i \Phi = \Phi_{i+1/2} - \Phi_{i-1/2}`
|
| :math:`~^{-i}` :
  :math:`\overline{\Phi}^i = ( \Phi_{i+1/2} + \Phi_{i-1/2} ) / 2`
|
| :math:`\delta_x` :
  :math:`\delta_x \Phi = \frac{1}{\Delta x} \delta_i \Phi`
|
| :math:`\overline{\nabla}` = horizontal gradient operator :
  :math:`\overline{\nabla} \Phi = \{ \delta_x \Phi , \delta_y \Phi \}`
|
| :math:`\overline{\nabla} \cdot` = horizontal divergence operator :
  :math:`\overline{\nabla}\cdot \vec{\mathrm{f}}  = 
  \frac{1}{\cal A} \{ \delta_i \Delta y \, \mathrm{f}_x 
                    + \delta_j \Delta x \, \mathrm{f}_y \}`
|
| :math:`\overline{\nabla}^2` = horizontal Laplacian operator :
  :math:`\overline{\nabla}^2 \Phi = 
     \overline{\nabla}\cdot \overline{\nabla}\Phi`
|

Time-stepping
=============

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
   time-stepping; 

#. as 1 but with an implicit linear free-surface;

#. as 1 or 2 but with variables staggered in time;

#. as 1 or 2 but with non-hydrostatic terms included;

#. as 2 or 3 but with non-linear free-surface.

In all the above configurations it is also possible to substitute the
Adams-Bashforth with an alternative time-stepping scheme for terms
evaluated explicitly in time. Since the over-arching algorithm is
independent of the particular time-stepping scheme chosen we will
describe first the over-arching algorithm, known as the pressure method,
with a rigid-lid model in :numref:`press_meth_rigid`. This
algorithm is essentially unchanged, apart for some coefficients, when
the rigid lid assumption is replaced with a linearized implicit
free-surface, described in :numref:`press_meth_linear`. These two flavors of the
pressure-method encompass all formulations of the model as it exists
today. The integration of explicit in time terms is out-lined in
:numref:`adams-bashforth` and put into the context of the overall algorithm
in :numref:`adams-bashforth-sync` and :numref:`adams-bashforth-staggered`.
Inclusion of non-hydrostatic terms
requires applying the pressure method in three dimensions instead of two
and this algorithm modification is described in
:numref:`non-hydrostatic`. Finally, the free-surface equation may be treated
more exactly, including non-linear terms, and this is described in
:numref:`nonlinear-freesurface`.

.. _press_meth_rigid:
 
Pressure method with rigid-lid
==============================

The horizontal momentum and continuity equations for the ocean
(:eq:`eq-ocean-mom` and :eq:`eq-ocean-cont`), or for the atmosphere
(:eq:`atmos-mom` and :eq:`atmos-cont`), can be summarized by:

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
   :label: rigid-lid-continuity

Here, :math:`H\widehat{u} = \int_H u dz` is the depth integral of
:math:`u`, similarly for :math:`H\widehat{v}`. The rigid-lid
approximation sets :math:`w=0` at the lid so that it does not move but
allows a pressure to be exerted on the fluid by the lid. The horizontal
momentum equations and vertically integrated continuity equation are be
discretized in time and space as follows:

.. math::
   u^{n+1} + \Delta t g \partial_x \eta^{n+1}
   =  u^{n} + \Delta t G_u^{(n+1/2)}
   :label: discrete-time-u

.. math::
   v^{n+1} + \Delta t g \partial_y \eta^{n+1}
   = v^{n} + \Delta t G_v^{(n+1/2)}
   :label: discrete-time-v

.. math::
   \partial_x H \widehat{u^{n+1}}
   + \partial_y H \widehat{v^{n+1}} = 0
   :label: discrete-time-cont-rigid-lid

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
:eq:`discrete-time-u`, :eq:`discrete-time-v` and
:eq:`discrete-time-cont-rigid-lid` can then be re-arranged as follows:

.. math:: u^{*} = u^{n} + \Delta t G_u^{(n+1/2)}
   :label: ustar-rigid-lid

.. math:: v^{*} = v^{n} + \Delta t G_v^{(n+1/2)}
   :label: vstar-rigid-lid

.. math:: \partial_x \Delta t g H \partial_x \eta^{n+1}
   + \partial_y \Delta t g H \partial_y \eta^{n+1}
   = \partial_x H \widehat{u^{*}}
   + \partial_y H \widehat{v^{*}} 
   :label: elliptic

.. math:: u^{n+1} = u^{*} - \Delta t g \partial_x \eta^{n+1}
   :label: un+1-rigid-lid

.. math:: v^{n+1} = v^{*} - \Delta t g \partial_y \eta^{n+1}
   :label: vn+1-rigid-lid

Equations :eq:`ustar-rigid-lid` to :eq:`vn+1-rigid-lid`, solved
sequentially, represent the pressure method algorithm used in the model.
The essence of the pressure method lies in the fact that any explicit
prediction for the flow would lead to a divergence flow field so a
pressure field must be found that keeps the flow non-divergent over each
step of the integration. The particular location in time of the pressure
field is somewhat ambiguous; in :numref:`pressure-method-rigid-lid` we
depicted as co-located with the future flow field (time level
:math:`n+1`) but it could equally have been drawn as staggered in time
with the flow.

  .. figure:: figs/pressure-method-rigid-lid.*
    :width: 70%
    :align: center
    :alt: pressure-method-rigid-lid
    :name: pressure-method-rigid-lid

    A schematic of the evolution in time of the pressure method algorithm. A prediction for the flow variables at time level :math:`n+1` is made based only on the explicit terms, :math:`G^{(n+^1/_2)}`, and denoted :math:`u^*`, :math:`v^*`. Next, a pressure field is found such that :math:`u^{n+1}`, :math:`v^{n+1}` will be non-divergent. Conceptually, the :math:`*` quantities exist at time level :math:`n+1` but they are intermediate and only temporary.


The correspondence to the code is as follows:

-  the prognostic phase, equations :eq:`ustar-rigid-lid` and
   :eq:`vstar-rigid-lid`, stepping forward :math:`u^n` and :math:`v^n` to
   :math:`u^{*}` and :math:`v^{*}` is coded in :filelink:`timestep.F <model/src/timestep.F>`

-  the vertical integration, :math:`H \widehat{u^*}` and :math:`H
   \widehat{v^*}`, divergence and inversion of the elliptic operator in
   equation :eq:`elliptic` is coded in :filelink:`solve_for_pressure.F <model/src/solve_for_pressure.F>`

-  finally, the new flow field at time level :math:`n+1` given by
   equations :eq:`un+1-rigid-lid` and :eq:`vn+1-rigid-lid` is calculated
   in :filelink:`correction_step.F <model/src/correction_step.F>`

The calling tree for these routines is as follows:

.. _call-tree-press-meth:

.. admonition:: Pressure method calling tree
  :class: note

    | :filelink:`FORWARD\_STEP <model/src/forward_step.F>`
    | :math:`\phantom{W}` :filelink:`DYNAMICS <model/src/dynamics.F>`
    | :math:`\phantom{WW}` :filelink:`TIMESTEP <model/src/timestep.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxxx}` :math:`u^*,v^*` :eq:`ustar-rigid-lid` , :eq:`vstar-rigid-lid`
    | :math:`\phantom{W}` :filelink:`SOLVE\_FOR\_PRESSURE <model/src/solve_for_pressure.F>`
    | :math:`\phantom{WW}` :filelink:`CALC\_DIV\_GHAT <model/src/calc_div_ghat.F>` :math:`\phantom{xxxxxxxxxxxxxxxx}` :math:`H\widehat{u^*},H\widehat{v^*}` :eq:`elliptic`
    | :math:`\phantom{WW}` :filelink:`CG2D <model/src/cg2d.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxxxxxxx}` :math:`\eta^{n+1}` :eq:`elliptic`
    | :math:`\phantom{W}` :filelink:`MOMENTUM\_CORRECTION\_STEP <model/src/momentum_correction_step.F>`
    | :math:`\phantom{WW}` :filelink:`CALC\_GRAD\_PHI\_SURF <model/src/calc_grad_phi_surf.F>` :math:`\phantom{xxxxxxxxxx}` :math:`\nabla \eta^{n+1}`
    | :math:`\phantom{WW}` :filelink:`CORRECTION\_STEP  <model/src/correction_step.F>` :math:`\phantom{xxxxxxxxxxxxw}` :math:`u^{n+1},v^{n+1}` :eq:`un+1-rigid-lid` , :eq:`vn+1-rigid-lid`

In general, the horizontal momentum time-stepping can contain some terms
that are treated implicitly in time, such as the vertical viscosity when
using the backward time-stepping scheme (:varlink:`implicitViscosity` =.TRUE.). The method used to solve
those implicit terms is provided in :numref:`implicit-backward-stepping`, and modifies equations
:eq:`discrete-time-u` and :eq:`discrete-time-v` to give:

.. math::

   \begin{aligned}
   u^{n+1} - \Delta t \partial_z A_v \partial_z u^{n+1}
   + \Delta t g \partial_x \eta^{n+1} & = & u^{n} + \Delta t G_u^{(n+1/2)}
   \\
   v^{n+1} - \Delta t \partial_z A_v \partial_z v^{n+1}
   + \Delta t g \partial_y \eta^{n+1} & = & v^{n} + \Delta t G_v^{(n+1/2)}\end{aligned}


.. _press_meth_linear:

Pressure method with implicit linear free-surface
=================================================

The rigid-lid approximation filters out external gravity waves
subsequently modifying the dispersion relation of barotropic Rossby
waves. The discrete form of the elliptic equation has some zero
eigen-values which makes it a potentially tricky or inefficient problem
to solve.

The rigid-lid approximation can be easily replaced by a linearization of
the free-surface equation which can be written:

.. math::
   \partial_t \eta + \partial_x H \widehat{u} + \partial_y H \widehat{v} = P-E+R
   :label: linear-free-surface=P-E

which differs from the depth integrated continuity equation with
rigid-lid (:eq:`rigid-lid-continuity`) by the time-dependent term and
fresh-water source term.

Equation :eq:`discrete-time-cont-rigid-lid` in the rigid-lid pressure
method is then replaced by the time discretization of
:eq:`linear-free-surface=P-E` which is:

.. math::
   \eta^{n+1}
   + \Delta t \partial_x H \widehat{u^{n+1}}
   + \Delta t \partial_y H \widehat{v^{n+1}}
   = \eta^{n} + \Delta t ( P - E )
   :label: discrete-time-backward-free-surface

where the use of flow at time level :math:`n+1` makes the method
implicit and backward in time. This is the preferred scheme since it
still filters the fast unresolved wave motions by damping them. A
centered scheme, such as Crank-Nicholson (see
:numref:`crank-nicolson_baro`), would alias the energy of the fast modes onto
slower modes of motion.

As for the rigid-lid pressure method, equations :eq:`discrete-time-u`,
:eq:`discrete-time-v` and :eq:`discrete-time-backward-free-surface` can be
re-arranged as follows:

.. math::
   u^{*} = u^{n} + \Delta t G_u^{(n+1/2)}
   :label: ustar-backward-free-surface

.. math::
   v^{*} = v^{n} + \Delta t G_v^{(n+1/2)}
   :label: vstar-backward-free-surface

.. math::
   \eta^* = \epsilon_{fs} ( \eta^{n} + \Delta t (P-E) )
            - \Delta t ( \partial_x H \widehat{u^{*}}
                            + \partial_y H \widehat{v^{*}} )
   :label: etastar-backward-free-surface

.. math::
   \partial_x g H \partial_x \eta^{n+1}
   + \partial_y g H \partial_y \eta^{n+1}
   - \frac{\epsilon_{fs} \eta^{n+1}}{\Delta t^2} =
   - \frac{\eta^*}{\Delta t^2}
   :label: elliptic-backward-free-surface

.. math::
   u^{n+1} = u^{*} - \Delta t g \partial_x \eta^{n+1}
   :label: un+1-backward-free-surface

.. math::
   v^{n+1} = v^{*} - \Delta t g \partial_y \eta^{n+1}
   :label: vn+1-backward-free-surface

Equations :eq:`ustar-backward-free-surface`
to :eq:`vn+1-backward-free-surface`, solved sequentially, represent the
pressure method algorithm with a backward implicit, linearized free
surface. The method is still formerly a pressure method because in the
limit of large :math:`\Delta t` the rigid-lid method is recovered.
However, the implicit treatment of the free-surface allows the flow to
be divergent and for the surface pressure/elevation to respond on a
finite time-scale (as opposed to instantly). To recover the rigid-lid
formulation, we introduced a switch-like parameter,
:math:`\epsilon_{fs}` (:varlink:`freesurfFac`), which selects between the free-surface and
rigid-lid; :math:`\epsilon_{fs}=1` allows the free-surface to evolve;
:math:`\epsilon_{fs}=0` imposes the rigid-lid. The evolution in time and
location of variables is exactly as it was for the rigid-lid model so
that :numref:`pressure-method-rigid-lid` is still applicable.
Similarly, the calling sequence, given :ref:`here <call-tree-press-meth>`, is as for the pressure-method.

.. _adams-bashforth:

Explicit time-stepping: Adams-Bashforth
=======================================

In describing the the pressure method above we deferred describing the
time discretization of the explicit terms. We have historically used the
quasi-second order Adams-Bashforth method for all explicit terms in both
the momentum and tracer equations. This is still the default mode of
operation but it is now possible to use alternate schemes for tracers
(see :numref:`tracer-advection`). In the previous sections, we summarized an explicit scheme as:

.. math::
   \tau^{*} = \tau^{n} + \Delta t G_\tau^{(n+1/2)}
   :label: taustar

where :math:`\tau` could be any prognostic variable (:math:`u`,
:math:`v`, :math:`\theta` or :math:`S`) and :math:`\tau^*` is an
explicit estimate of :math:`\tau^{n+1}` and would be exact if not for
implicit-in-time terms. The parenthesis about :math:`n+1/2` indicates
that the term is explicit and extrapolated forward in time and for this
we use the quasi-second order Adams-Bashforth method:

.. math::
   G_\tau^{(n+1/2)} = ( 3/2 + \epsilon_{AB}) G_\tau^n
   - ( 1/2 + \epsilon_{AB}) G_\tau^{n-1}
   :label: adams-bashforth2

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
terms outside the Adams-Bashforth extrapolation, by turning off the logical flag :varlink:`forcing_In_AB`
(parameter file ``data``, namelist ``PARM01``, default value = TRUE) and then setting :varlink:`tracForcingOutAB`
(default=0), :varlink:`momForcingOutAB` (default=0), and :varlink:`momDissip_In_AB` (parameter file ``data``, namelist ``PARM01``,
default value = TRUE), respectively for the tracer terms, momentum forcing terms, and the dissipation terms.

A stability analysis for an oscillation equation should be given at this
point.

A stability analysis for a relaxation equation should be given at this
point.


  .. figure:: figs/oscil+damp_AB2.*
    :width: 100%
    :align: center
    :alt: stability_analysis
    :name: oscil+damp_AB2

    Oscillatory and damping response of quasi-second order Adams-Bashforth scheme for different values of the  :math:`\epsilon _{AB}` parameter (0.0, 0.1, 0.25, from top to bottom) The analytical solution (in black), the physical mode (in blue) and the numerical mode (in red) are represented with a CFL step of 0.1. The left column represents the oscillatory response on the complex plane for CFL ranging from 0.1 up to 0.9. The right column represents the damping response amplitude (y-axis) function of the CFL (x-axis).

.. _implicit-backward-stepping:

Implicit time-stepping: backward method
=======================================

Vertical diffusion and viscosity can be treated implicitly in time using
the backward method which is an intrinsic scheme. Recently, the option
to treat the vertical advection implicitly has been added, but not yet
tested; therefore, the description hereafter is limited to diffusion and
viscosity. For tracers, the time discretized equation is:

.. math::
   \tau^{n+1} - \Delta t \partial_r \kappa_v \partial_r \tau^{n+1} = \tau^{n} + \Delta t G_\tau^{(n+1/2)}
   :label: implicit-diffusion

where :math:`G_\tau^{(n+1/2)}` is the remaining explicit terms
extrapolated using the Adams-Bashforth method as described above.
Equation :eq:`implicit-diffusion` can be split split into:

.. math::
   \tau^* = \tau^{n} + \Delta t G_\tau^{(n+1/2)}
   :label: taustar-implicit

.. math::
   \tau^{n+1} = {\cal L}_\tau^{-1} ( \tau^* )
   :label: tau-n+1-implicit

where :math:`{\cal L}_\tau^{-1}` is the inverse of the operator

.. math:: {\cal L}_\tau = \left[ 1 + \Delta t \partial_r \kappa_v \partial_r \right]

Equation :eq:`taustar-implicit` looks exactly as :eq:`taustar` while
:eq:`tau-n+1-implicit` involves an operator or matrix inversion. By
re-arranging :eq:`implicit-diffusion` in this way we have cast the method
as an explicit prediction step and an implicit step allowing the latter
to be inserted into the over all algorithm with minimal interference.

The calling sequence for stepping forward a tracer variable such as temperature with implicit diffusion is
as follows:

.. _adams-bash-calltree:

.. admonition:: Adams-Bashforth calling tree
  :class: note

    | :filelink:`FORWARD\_STEP <model/src/forward_step.F>`
    | :math:`\phantom{W}` :filelink:`THERMODYNAMICS <model/src/thermodynamics.F>`
    | :math:`\phantom{WW}` :filelink:`TEMP\_INTEGRATE <model/src/temp_integrate.F>`
    | :math:`\phantom{WWW}` :filelink:`GAD\_CALC\_RHS <pkg/generic_advdiff/gad_calc_rhs.F>` :math:`\phantom{xxxxxxxxxw}` :math:`G_\theta^n = G_\theta( u, \theta^n)`
    | :math:`\phantom{WWW}` either
    | :math:`\phantom{WWWW}` :filelink:`EXTERNAL\_FORCING <model/src/external_forcing.F>` :math:`\phantom{xxxx}` :math:`G_\theta^n = G_\theta^n + {\cal Q}`
    | :math:`\phantom{WWWW}` :filelink:`ADAMS\_BASHFORTH2 <model/src/adams_bashforth2.F>` :math:`\phantom{xxi}` :math:`G_\theta^{(n+1/2)}` :eq:`adams-bashforth2`
    | :math:`\phantom{WWW}` or
    | :math:`\phantom{WWWW}` :filelink:`EXTERNAL\_FORCING <model/src/external_forcing.F>` :math:`\phantom{xxxx}` :math:`G_\theta^{(n+1/2)} = G_\theta^{(n+1/2)} + {\cal Q}`
    | :math:`\phantom{WW}` :filelink:`TIMESTEP\_TRACER <model/src/timestep_tracer.F>` :math:`\phantom{xxxxxxxxxx}` :math:`\tau^*` :eq:`taustar`
    | :math:`\phantom{WW}` :filelink:`IMPLDIFF  <model/src/impldiff.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxw}` :math:`\tau^{(n+1)}` :eq:`tau-n+1-implicit`

In order to fit within the pressure method, the implicit viscosity must
not alter the barotropic flow. In other words, it can only redistribute
momentum in the vertical. The upshot of this is that although vertical
viscosity may be backward implicit and unconditionally stable, no-slip
boundary conditions may not be made implicit and are thus cast as a an
explicit drag term.

.. _adams-bashforth-sync:

Synchronous time-stepping: variables co-located in time
=======================================================

  .. figure:: figs/adams-bashforth-sync.*
    :width: 70%
    :align: center
    :alt: adams-bash-sync
    :name: adams-bash-sync

    A schematic of the explicit Adams-Bashforth and implicit time-stepping phases of the algorithm. All prognostic variables are co-located in time. Explicit tendencies are evaluated at time level :math:`n` as a function of the state at that time level (dotted arrow). The explicit tendency from the previous time level, :math:`n-1`, is used to extrapolate tendencies to :math:`n+1/2` (dashed arrow). This extrapolated tendency allows variables to be stably integrated forward-in-time to render an estimate (:math:`*` -variables) at the :math:`n+1` time level (solid arc-arrow). The operator :math:`{\cal L}` formed from implicit-in-time terms is solved to yield the state variables at time level :math:`n+1`.


The Adams-Bashforth extrapolation of explicit tendencies fits neatly
into the pressure method algorithm when all state variables are
co-located in time. The algorithm can be represented by the sequential solution of the
follow equations:

.. math::   
   G_{\theta,S}^{n} = G_{\theta,S} ( u^n, \theta^n, S^n )
   :label: Gt-n-sync

.. math::
   G_{\theta,S}^{(n+1/2)} = (3/2+\epsilon_{AB}) G_{\theta,S}^{n}-(1/2+\epsilon_{AB}) G_{\theta,S}^{n-1}
   :label: Gt-n+5-sync

.. math::
   (\theta^*,S^*) = (\theta^{n},S^{n}) + \Delta t G_{\theta,S}^{(n+1/2)}
   :label: tstar-sync

.. math::
   (\theta^{n+1},S^{n+1}) = {\cal L}^{-1}_{\theta,S} (\theta^*,S^*)
   :label: t-n+1-sync

.. math::
   \phi^n_{hyd} = \int b(\theta^n,S^n) dr
   :label: phi-hyd-sync

.. math::
   \vec{\bf G}_{\vec{\bf v}}^{n} = \vec{\bf G}_{\vec{\bf v}} ( \vec{\bf v}^n, \phi^n_{hyd} )
   :label: Gv-n-sync

.. math::
   \vec{\bf G}_{\vec{\bf v}}^{(n+1/2)} = (3/2 + \epsilon_{AB} ) \vec{\bf G}_{\vec{\bf v}}^{n} - (1/2 + \epsilon_{AB} ) \vec{\bf G}_{\vec{\bf v}}^{n-1}
   :label: Gv-n+5-sync

.. math::
   \vec{\bf v}^{*} = \vec{\bf v}^{n} + \Delta t \vec{\bf G}_{\vec{\bf v}}^{(n+1/2)}
   :label: vstar-sync

.. math::
   \vec{\bf v}^{**} = {\cal L}_{\vec{\bf v}}^{-1} ( \vec{\bf v}^* )
   :label: vstarstar-sync

.. math::
   \eta^* = \epsilon_{fs} \left( \eta^{n} + \Delta t (P-E) \right)- \Delta t
     \nabla \cdot H \widehat{ \vec{\bf v}^{**} }
   :label: nstar-sync

.. math::
   \nabla \cdot g H \nabla \eta^{n+1} - \frac{\epsilon_{fs} \eta^{n+1}}{\Delta t^2} ~ = ~ - \frac{\eta^*}{\Delta t^2}
   :label: elliptic-sync

.. math::
   \vec{\bf v}^{n+1} = \vec{\bf v}^{**} - \Delta t g \nabla \eta^{n+1}
   :label: v-n+1-sync

:numref:`adams-bash-sync` illustrates the location of variables
in time and evolution of the algorithm with time. The Adams-Bashforth
extrapolation of the tracer tendencies is illustrated by the dashed
arrow, the prediction at :math:`n+1` is indicated by the solid arc.
Inversion of the implicit terms, :math:`{\cal
L}^{-1}_{\theta,S}`, then yields the new tracer fields at :math:`n+1`.
All these operations are carried out in subroutine :filelink:`THERMODYNAMICS <model/src/thermodynamics.F>` and
subsidiaries, which correspond to equations :eq:`Gt-n-sync` to
:eq:`t-n+1-sync`. Similarly illustrated is the Adams-Bashforth
extrapolation of accelerations, stepping forward and solving of implicit
viscosity and surface pressure gradient terms, corresponding to
equations :eq:`Gv-n-sync` to :eq:`v-n+1-sync`. These operations are
carried out in subroutines :filelink:`DYNAMICS <model/src/dynamics.F>`,
:filelink:`SOLVE\_FOR\_PRESSURE <model/src/solve_for_pressure.F>` and
:filelink:`MOMENTUM\_CORRECTION\_STEP <model/src/momentum_correction_step.F>`.
This, then, represents an entire algorithm
for stepping forward the model one time-step. The corresponding calling
tree for the overall synchronous algorithm using
Adams-Bashforth time-stepping is given below. The place where the model geometry
**hFac** factors) is updated is added here but is only relevant
for the non-linear free-surface algorithm.
For completeness, the external forcing,
ocean and atmospheric physics have been added, although they are mainly optional.

.. admonition:: Synchronous Adams-Bashforth calling tree
  :class: note

    | :filelink:`FORWARD\_STEP <model/src/forward_step.F>`
    | :math:`\phantom{WWW}` :filelink:`EXTERNAL\_FIELDS\_LOAD <model/src/external_fields_load.F>`
    | :math:`\phantom{WWW}` :filelink:`DO\_ATMOSPHERIC\_PHYS <model/src/do_atmospheric_phys.F>`
    | :math:`\phantom{WWW}` :filelink:`DO\_OCEANIC\_PHYS <model/src/do_oceanic_phys.F>`
    | :math:`\phantom{WW}` :filelink:`THERMODYNAMICS <model/src/thermodynamics.F>`
    | :math:`\phantom{WWW}` :filelink:`CALC\_GT <model/src/calc_gt.F>`
    | :math:`\phantom{WWWW}` :filelink:`GAD\_CALC\_RHS <pkg/generic_advdiff/gad_calc_rhs.F>` :math:`\phantom{xxxxxxxxxxxxxlwww}` :math:`G_\theta^n = G_\theta( u, \theta^n )` :eq:`Gt-n-sync`
    | :math:`\phantom{WWWW}` :filelink:`EXTERNAL\_FORCING <model/src/external_forcing.F>` :math:`\phantom{xxxxxxxxxxlww}` :math:`G_\theta^n = G_\theta^n + {\cal Q}`
    | :math:`\phantom{WWWW}` :filelink:`ADAMS\_BASHFORTH2 <model/src/adams_bashforth2.F>` :math:`\phantom{xxxxxxxxxxxw}` :math:`G_\theta^{(n+1/2)}` :eq:`Gt-n+5-sync`
    | :math:`\phantom{WWW}` :filelink:`TIMESTEP\_TRACER <model/src/timestep_tracer.F>` :math:`\phantom{xxxxxxxxxxxxxxxww}` :math:`\theta^*` :eq:`tstar-sync`
    | :math:`\phantom{WWW}` :filelink:`IMPLDIFF  <model/src/impldiff.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxxvwww}` :math:`\theta^{(n+1)}` :eq:`t-n+1-sync`
    | :math:`\phantom{WW}` :filelink:`DYNAMICS <model/src/dynamics.F>`
    | :math:`\phantom{WWW}` :filelink:`CALC\_PHI\_HYD <model/src/calc_phi_hyd.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxi}` :math:`\phi_{hyd}^n` :eq:`phi-hyd-sync`
    | :math:`\phantom{WWW}` :filelink:`MOM\_FLUXFORM <pkg/mom_fluxform/mom_fluxform.F>` or :filelink:`MOM\_VECINV <pkg/mom_vecinv/mom_vecinv.F>` :math:`\phantom{xxi}` :math:`G_{\vec{\bf v}}^n` :eq:`Gv-n-sync`
    | :math:`\phantom{WWW}` :filelink:`TIMESTEP <model/src/timestep.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxxxxxxx}` :math:`\vec{\bf v}^*` :eq:`Gv-n+5-sync`, :eq:`vstar-sync`
    | :math:`\phantom{WWW}` :filelink:`IMPLDIFF  <model/src/impldiff.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxxxxxlw}` :math:`\vec{\bf v}^{**}` :eq:`vstarstar-sync`
    | :math:`\phantom{WW}` :filelink:`UPDATE\_R\_STAR <model/src/update_r_star.F>` or :filelink:`UPDATE\_SURF\_DR <model/src/update_surf_dr.F>` (NonLin-FS only)
    | :math:`\phantom{WW}` :filelink:`SOLVE\_FOR\_PRESSURE <model/src/solve_for_pressure.F>`
    | :math:`\phantom{WWW}` :filelink:`CALC\_DIV\_GHAT <model/src/calc_div_ghat.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxx}` :math:`\eta^*` :eq:`nstar-sync`
    | :math:`\phantom{WWW}` :filelink:`CG2D <model/src/cg2d.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxxxxxxxxxxi}` :math:`\eta^{n+1}` :eq:`elliptic-sync`
    | :math:`\phantom{WW}` :filelink:`MOMENTUM\_CORRECTION\_STEP <model/src/momentum_correction_step.F>`
    | :math:`\phantom{WWW}` :filelink:`CALC\_GRAD\_PHI\_SURF <model/src/calc_grad_phi_surf.F>` :math:`\phantom{xxxxxxxxxxxxxx}` :math:`\nabla \eta^{n+1}`
    | :math:`\phantom{WWW}` :filelink:`CORRECTION\_STEP  <model/src/correction_step.F>` :math:`\phantom{xxxxxxxxxxxxxxxxw}` :math:`u^{n+1},v^{n+1}` :eq:`v-n+1-sync`
    | :math:`\phantom{WW}` :filelink:`TRACERS\_CORRECTION\_STEP <model/src/tracers_correction_step.F>`
    | :math:`\phantom{WWW}` :filelink:`CYCLE\_TRACER <model/src/cycle_tracer.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxx}` :math:`\theta^{n+1}`
    | :math:`\phantom{WWW}` :filelink:`SHAP_FILT_APPLY_TS  <pkg/shap_filt/shap_filt_apply_ts.F>` or :filelink:`ZONAL_FILT_APPLY_TS  <pkg/zonal_filt/zonal_filt_apply_ts.F>`
    | :math:`\phantom{WWW}` :filelink:`CONVECTIVE_ADJUSTMENT  <model/src/convective_adjustment.F>`


.. _adams-bashforth-staggered:

Staggered baroclinic time-stepping
==================================

  .. figure:: figs/adams-bashforth-staggered.*
    :width: 80%
    :align: center
    :alt: adams-bash-staggered
    :name: adams-bash-staggered

    A schematic of the explicit Adams-Bashforth and implicit time-stepping phases of the algorithm but with staggering in time of thermodynamic variables with the flow. Explicit momentum tendencies are evaluated at time level :math:`n-1/2` as a function of the flow field at that time level :math:`n-1/2`. The explicit tendency from the previous time level, :math:`n-3/2`, is used to extrapolate tendencies to :math:`n` (dashed arrow). The hydrostatic pressure/geo-potential  :math:`\phi _{hyd}` is evaluated directly at time level :math:`n` (vertical arrows) and used with the extrapolated tendencies to step forward the flow variables from :math:`n-1/2` to :math:`n+1/2` (solid arc-arrow). The implicit-in-time operator  :math:`{\cal L}_{\bf u,v}` (vertical arrows) is then applied to the previous estimation of the the flow field (:math:`*` -variables) and yields to the two velocity components :math:`u,v` at time level :math:`n+1/2`. These are then used to calculate the advection term (dashed arc-arrow) of the thermo-dynamics tendencies at time step :math:`n`. The extrapolated thermodynamics tendency, from time level :math:`n-1` and :math:`n` to :math:`n+1/2`, allows thermodynamic variables to be stably integrated forward-in-time (solid arc-arrow) up to time level :math:`n+1`.

For well-stratified problems, internal gravity waves may be the limiting
process for determining a stable time-step. In the circumstance, it is
more efficient to stagger in time the thermodynamic variables with the
flow variables. :numref:`adams-bash-staggered` illustrates the
staggering and algorithm. The key difference between this and
:numref:`adams-bash-sync` is that the thermodynamic variables are
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
entire algorithm, :eq:`Gt-n-sync` to :eq:`v-n+1-sync`, annotating the
position in time of variables appropriately:

.. math::
   \phi^{n}_{hyd} =  \int b(\theta^{n},S^{n}) dr
   :label: phi-hyd-staggered

.. math::
   \vec{\bf G}_{\vec{\bf v}}^{n-1/2}  =  \vec{\bf G}_{\vec{\bf v}} ( \vec{\bf v}^{n-1/2} )
   :label: Gv-n-staggered

.. math::
   \vec{\bf G}_{\vec{\bf v}}^{(n)} =  (3/2 + \epsilon_{AB} ) \vec{\bf G}_{\vec{\bf v}}^{n-1/2} - (1/2 + \epsilon_{AB} ) \vec{\bf G}_{\vec{\bf v}}^{n-3/2}
   :label: Gv-n+5-staggered

.. math::
   \vec{\bf v}^{*}  =  \vec{\bf v}^{n-1/2} + \Delta t \left( \vec{\bf G}_{\vec{\bf v}}^{(n)} - \nabla \phi_{hyd}^{n} \right)
   :label: vstar-staggered

.. math::
   \vec{\bf v}^{**}  =  {\cal L}_{\vec{\bf v}}^{-1} ( \vec{\bf v}^* )
   :label: vstarstar-staggered

.. math::
   \eta^*  = \epsilon_{fs} \left( \eta^{n-1/2} + \Delta t (P-E)^n \right)- \Delta t
     \nabla \cdot H \widehat{ \vec{\bf v}^{**} }
   :label: nstar-staggered

.. math::
   \nabla \cdot g H \nabla \eta^{n+1/2}  -  \frac{\epsilon_{fs} \eta^{n+1/2}}{\Delta t^2}
   ~ = ~ - \frac{\eta^*}{\Delta t^2}
   :label: elliptic-staggered

.. math::
   \vec{\bf v}^{n+1/2}  =  \vec{\bf v}^{**} - \Delta t g \nabla \eta^{n+1/2}
   :label: v-n+1-staggered

.. math::
   G_{\theta,S}^{n}  =  G_{\theta,S} ( u^{n+1/2}, \theta^{n}, S^{n} )
   :label: Gt-n-staggered

.. math::
   G_{\theta,S}^{(n+1/2)}  =  (3/2+\epsilon_{AB}) G_{\theta,S}^{n}-(1/2+\epsilon_{AB}) G_{\theta,S}^{n-1}
   :label: Gt-n+5-staggered

.. math::
   (\theta^*,S^*)  =  (\theta^{n},S^{n}) + \Delta t G_{\theta,S}^{(n+1/2)}
   :label: tstar-staggered

.. math::
   (\theta^{n+1},S^{n+1})  =  {\cal L}^{-1}_{\theta,S} (\theta^*,S^*)
   :label: t-n+1-staggered

The corresponding calling tree is given below. The staggered algorithm is
activated with the run-time flag :varlink:`staggerTimeStep` =.TRUE. in
parameter file ``data``, namelist ``PARM01``.

.. admonition:: Staggered Adams-Bashforth calling tree
  :class: note

    | :filelink:`FORWARD\_STEP <model/src/forward_step.F>`
    | :math:`\phantom{WWW}` :filelink:`EXTERNAL\_FIELDS\_LOAD <model/src/external_fields_load.F>`
    | :math:`\phantom{WWW}` :filelink:`DO\_ATMOSPHERIC\_PHYS <model/src/do_atmospheric_phys.F>`
    | :math:`\phantom{WWW}` :filelink:`DO\_OCEANIC\_PHYS <model/src/do_oceanic_phys.F>`
    | :math:`\phantom{WW}` :filelink:`DYNAMICS <model/src/dynamics.F>`
    | :math:`\phantom{WWW}` :filelink:`CALC\_PHI\_HYD <model/src/calc_phi_hyd.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxi}` :math:`\phi_{hyd}^n` :eq:`phi-hyd-staggered`
    | :math:`\phantom{WWW}` :filelink:`MOM\_FLUXFORM <pkg/mom_fluxform/mom_fluxform.F>` or :filelink:`MOM\_VECINV <pkg/mom_vecinv/mom_vecinv.F>` :math:`\phantom{xxi}` :math:`G_{\vec{\bf v}}^{n-1/2}` :eq:`Gv-n-staggered`
    | :math:`\phantom{WWW}` :filelink:`TIMESTEP <model/src/timestep.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxxxxxxx}` :math:`\vec{\bf v}^*` :eq:`Gv-n+5-staggered`, :eq:`vstar-staggered`
    | :math:`\phantom{WWW}` :filelink:`IMPLDIFF  <model/src/impldiff.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxxxxxlw}` :math:`\vec{\bf v}^{**}` :eq:`vstarstar-staggered`
    | :math:`\phantom{WW}` :filelink:`UPDATE\_R\_STAR <model/src/update_r_star.F>` or :filelink:`UPDATE\_SURF\_DR <model/src/update_surf_dr.F>` (NonLin-FS only)
    | :math:`\phantom{WW}` :filelink:`SOLVE\_FOR\_PRESSURE <model/src/solve_for_pressure.F>`
    | :math:`\phantom{WWW}` :filelink:`CALC\_DIV\_GHAT <model/src/calc_div_ghat.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxx}` :math:`\eta^*` :eq:`nstar-staggered`
    | :math:`\phantom{WWW}` :filelink:`CG2D <model/src/cg2d.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxxxxxxxxxxi}` :math:`\eta^{n+1/2}` :eq:`elliptic-staggered`
    | :math:`\phantom{WW}` :filelink:`MOMENTUM\_CORRECTION\_STEP <model/src/momentum_correction_step.F>`
    | :math:`\phantom{WWW}` :filelink:`CALC\_GRAD\_PHI\_SURF <model/src/calc_grad_phi_surf.F>` :math:`\phantom{xxxxxxxxxxxxxx}` :math:`\nabla \eta^{n+1/2}`
    | :math:`\phantom{WWW}` :filelink:`CORRECTION\_STEP  <model/src/correction_step.F>` :math:`\phantom{xxxxxxxxxxxxxxxxw}` :math:`u^{n+1/2},v^{n+1/2}` :eq:`v-n+1-staggered`
    | :math:`\phantom{WW}` :filelink:`THERMODYNAMICS <model/src/thermodynamics.F>`
    | :math:`\phantom{WWW}` :filelink:`CALC\_GT <model/src/calc_gt.F>`
    | :math:`\phantom{WWWW}` :filelink:`GAD\_CALC\_RHS <pkg/generic_advdiff/gad_calc_rhs.F>` :math:`\phantom{xxxxxxxxxxxxxlwww}` :math:`G_\theta^n = G_\theta( u, \theta^n )` :eq:`Gt-n-staggered`
    | :math:`\phantom{WWWW}` :filelink:`EXTERNAL\_FORCING <model/src/external_forcing.F>` :math:`\phantom{xxxxxxxxxxlww}` :math:`G_\theta^n = G_\theta^n + {\cal Q}`
    | :math:`\phantom{WWWW}` :filelink:`ADAMS\_BASHFORTH2 <model/src/adams_bashforth2.F>` :math:`\phantom{xxxxxxxxxxxw}` :math:`G_\theta^{(n+1/2)}` :eq:`Gt-n+5-staggered`
    | :math:`\phantom{WWW}` :filelink:`TIMESTEP\_TRACER <model/src/timestep_tracer.F>` :math:`\phantom{xxxxxxxxxxxxxxxww}` :math:`\theta^*` :eq:`tstar-staggered`
    | :math:`\phantom{WWW}` :filelink:`IMPLDIFF  <model/src/impldiff.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxxvwww}` :math:`\theta^{(n+1)}` :eq:`t-n+1-staggered`
    | :math:`\phantom{WW}` :filelink:`TRACERS\_CORRECTION\_STEP <model/src/tracers_correction_step.F>`
    | :math:`\phantom{WWW}` :filelink:`CYCLE\_TRACER <model/src/cycle_tracer.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxx}` :math:`\theta^{n+1}`
    | :math:`\phantom{WWW}` :filelink:`SHAP_FILT_APPLY_TS  <pkg/shap_filt/shap_filt_apply_ts.F>` or :filelink:`ZONAL_FILT_APPLY_TS  <pkg/zonal_filt/zonal_filt_apply_ts.F>`
    | :math:`\phantom{WWW}` :filelink:`CONVECTIVE_ADJUSTMENT  <model/src/convective_adjustment.F>`

The only difficulty with this approach is apparent in equation
:eq:`Gt-n-staggered` and illustrated by the dotted arrow connecting
:math:`u,v^{n+1/2}` with :math:`G_\theta^{n}`. The flow used to advect
tracers around is not naturally located in time. This could be avoided
by applying the Adams-Bashforth extrapolation to the tracer field itself
and advecting that around but this approach is not yet available. We’re
not aware of any detrimental effect of this feature. The difficulty lies
mainly in interpretation of what time-level variables and terms
correspond to.

.. _non-hydrostatic:

Non-hydrostatic formulation
===========================


The non-hydrostatic formulation re-introduces the full vertical momentum
equation and requires the solution of a 3-D elliptic equations for
non-hydrostatic pressure perturbation. We still integrate vertically for
the hydrostatic pressure and solve a 2-D elliptic equation for the
surface pressure/elevation for this reduces the amount of work needed to
solve for the non-hydrostatic pressure.

The momentum equations are discretized in time as follows:

.. math::
   \frac{1}{\Delta t} u^{n+1} + g \partial_x \eta^{n+1} + \partial_x \phi_{nh}^{n+1}
   = \frac{1}{\Delta t} u^{n} + G_u^{(n+1/2)}
   :label: discrete-time-u-nh

.. math::
   \frac{1}{\Delta t} v^{n+1} + g \partial_y \eta^{n+1} + \partial_y \phi_{nh}^{n+1}
   = \frac{1}{\Delta t} v^{n} + G_v^{(n+1/2)}
   :label: discrete-time-v-nh

.. math::
   \frac{1}{\Delta t} w^{n+1} + \partial_r \phi_{nh}^{n+1}
   = \frac{1}{\Delta t} w^{n} + G_w^{(n+1/2)}
   :label: discrete-time-w-nh

which must satisfy the discrete-in-time depth integrated continuity,
equation :eq:`discrete-time-backward-free-surface` and the local
continuity equation

.. math::
   \partial_x u^{n+1} + \partial_y v^{n+1} + \partial_r w^{n+1} = 0
   :label: non-divergence-nh

As before, the explicit predictions for momentum are consolidated as:

.. math::

   \begin{aligned}
   u^* & = & u^n + \Delta t G_u^{(n+1/2)} \\
   v^* & = & v^n + \Delta t G_v^{(n+1/2)} \\
   w^* & = & w^n + \Delta t G_w^{(n+1/2)}\end{aligned}

but this time we introduce an intermediate step by splitting the
tendency of the flow as follows:

.. math::

   \begin{aligned}
   u^{n+1} = u^{**} - \Delta t \partial_x \phi_{nh}^{n+1}
   & &
   u^{**} = u^{*} - \Delta t g \partial_x \eta^{n+1} \\
   v^{n+1} = v^{**} - \Delta t \partial_y \phi_{nh}^{n+1}
   & &
   v^{**} = v^{*} - \Delta t g \partial_y \eta^{n+1}\end{aligned}

Substituting into the depth integrated continuity
(equation :eq:`discrete-time-backward-free-surface`) gives

.. math::
   \partial_x H \partial_x \left( g \eta^{n+1} + \widehat{\phi}_{nh}^{n+1} \right)
   + \partial_y H \partial_y \left( g \eta^{n+1} + \widehat{\phi}_{nh}^{n+1} \right)
    - \frac{\epsilon_{fs}\eta^{n+1}}{\Delta t^2}
   = - \frac{\eta^*}{\Delta t^2}
   :label: substituting-in-cont

which is approximated by equation :eq:`elliptic-backward-free-surface`
on the basis that i) :math:`\phi_{nh}^{n+1}` is not yet known and ii)
:math:`\nabla \widehat{\phi}_{nh}
<< g \nabla \eta`. If :eq:`elliptic-backward-free-surface` is solved
accurately then the implication is that :math:`\widehat{\phi}_{nh}
\approx 0` so that the non-hydrostatic pressure field does not drive
barotropic motion.

The flow must satisfy non-divergence (equation :eq:`non-divergence-nh`)
locally, as well as depth integrated, and this constraint is used to
form a 3-D elliptic equations for :math:`\phi_{nh}^{n+1}`:

.. math::
   \partial_{xx} \phi_{nh}^{n+1} + \partial_{yy} \phi_{nh}^{n+1} +
   \partial_{rr} \phi_{nh}^{n+1} =
   \partial_x u^{**} + \partial_y v^{**} + \partial_r w^{*}
   :label: elliptic-pnh

The entire algorithm can be summarized as the sequential solution of the
following equations:

.. math::
   u^{*} = u^{n} + \Delta t G_u^{(n+1/2)}
   :label: ustar-nh

.. math::   
   v^{*} = v^{n} + \Delta t G_v^{(n+1/2)}
   :label: vstar-nh
 
.. math::
   w^{*} = w^{n} + \Delta t G_w^{(n+1/2)}
   :label: wstar-nh

.. math::
   \eta^* ~ = ~ \epsilon_{fs} \left( \eta^{n} + \Delta t (P-E) \right)
   - \Delta t \left( \partial_x H \widehat{u^{*}}
                       + \partial_y H \widehat{v^{*}} \right)
   :label: etastar-nh

.. math::
    \partial_x g H \partial_x \eta^{n+1}
   + \partial_y g H \partial_y \eta^{n+1}
   - \frac{\epsilon_{fs} \eta^{n+1}}{\Delta t^2}
   ~ = ~ - \frac{\eta^*}{\Delta t^2}
   :label: elliptic-nh

.. math::
   u^{**} = u^{*} - \Delta t g \partial_x \eta^{n+1}
   :label: unx-nh

.. math::
   v^{**} = v^{*} - \Delta t g \partial_y \eta^{n+1}
   :label: vnx-nh

.. math::
   \partial_{xx} \phi_{nh}^{n+1} + \partial_{yy} \phi_{nh}^{n+1} +
   \partial_{rr} \phi_{nh}^{n+1} =
   \partial_x u^{**} + \partial_y v^{**} + \partial_r w^{*}
   :label: phi-nh

.. math::
   u^{n+1} = u^{**} - \Delta t \partial_x \phi_{nh}^{n+1}
   :label: un+1-nh

.. math::
   v^{n+1} = v^{**} - \Delta t \partial_y \phi_{nh}^{n+1}
   :label: vn+1-nh}

.. math::
   \partial_r w^{n+1} = - \partial_x u^{n+1} - \partial_y v^{n+1}
   :label: wn+1-nh

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
:eq:`discrete-time-backward-free-surface`, and assuming no
non-hydrostatic effects (:math:`\epsilon_{nh} = 0`) is:

.. math::
   \epsilon_{fs} {\eta}^{n+1} -
   {\bf \nabla}_h \cdot \Delta t^2 (R_o-R_{fixed}) {\bf \nabla}_h b_s
   {\eta}^{n+1} = {\eta}^*
   :label: eq-solve2D

where

.. math::
   {\eta}^* = \epsilon_{fs} \: {\eta}^{n} -
   \Delta t {\bf \nabla}_h \cdot \int_{R_{fixed}}^{R_o} \vec{\bf v}^* dr
   \: + \: \epsilon_{fw} \Delta t (P-E)^{n}
   :label: eq-solve2D_rhs

.. admonition:: S/R  :filelink:`SOLVE_FOR_PRESSURE <model/src/solve_for_pressure.F>`
  :class: note

    | :math:`u^*` : **gU** ( :filelink:`DYNVARS.h <model/inc/DYNVARS.h>` )
    | :math:`v^*` : **gV** ( :filelink:`DYNVARS.h <model/inc/DYNVARS.h>` )
    | :math:`{\eta}^*` : **cg2d_b** ( :filelink:`SOLVE_FOR_PRESSURE.h <model/inc/SOLVE_FOR_PRESSURE.h>` )
    | :math:`{\eta}^{n+1}` : **etaN** ( :filelink:`DYNVARS.h <model/inc/DYNVARS.h>` )


Once :math:`{\eta}^{n+1}` has been found, substituting into
:eq:`discrete-time-u`, :eq:`discrete-time-v` yields
:math:`\vec{\bf v}^{n+1}` if the model is hydrostatic
(:math:`\epsilon_{nh}=0`):

.. math::
   \vec{\bf v}^{n+1} = \vec{\bf v}^{*}
   - \Delta t {\bf \nabla}_h b_s {\eta}^{n+1}

This is known as the correction step. However, when the model is
non-hydrostatic (:math:`\epsilon_{nh}=1`) we need an additional step and
an additional equation for :math:`\phi'_{nh}`. This is obtained by
substituting :eq:`discrete-time-u-nh`, :eq:`discrete-time-v-nh` and
:eq:`discrete-time-w-nh` into continuity:

.. math::
   [ {\bf \nabla}_h^2 + \partial_{rr} ] {\phi'_{nh}}^{n+1}
   = \frac{1}{\Delta t} 
   {\bf \nabla}_h \cdot \vec{\bf v}^{**} + \partial_r \dot{r}^*
   :label: sub-u-v-w-in-cont 

where

.. math:: \vec{\bf v}^{**} = \vec{\bf v}^* - \Delta t {\bf \nabla}_h b_s {\eta}^{n+1}

Note that :math:`\eta^{n+1}` is also used to update the second RHS term
:math:`\partial_r \dot{r}^*` since the vertical velocity at the surface
(:math:`\dot{r}_{surf}`) is evaluated as
:math:`(\eta^{n+1} - \eta^n) / \Delta t`.

Finally, the horizontal velocities at the new time level are found by:

.. math::
   \vec{\bf v}^{n+1} = \vec{\bf v}^{**}
   - \epsilon_{nh} \Delta t {\bf \nabla}_h {\phi'_{nh}}^{n+1}
   :label: v-new-time-lev 

and the vertical velocity is found by integrating the continuity
equation vertically. Note that, for the convenience of the restart
procedure, the vertical integration of the continuity equation has been
moved to the beginning of the time step (instead of at the end), without
any consequence on the solution.

.. admonition:: S/R  :filelink:`CORRECTION_STEP <model/src/correction_step.F>`
  :class: note

    | :math:`{\eta}^{n+1}` : **etaN** ( :filelink:`DYNVARS.h <model/inc/DYNVARS.h>` )
    | :math:`{\phi}^{n+1}_{nh}` : **phi_nh** ( :filelink:`NH_VARS.h <model/inc/NH_VARS.h>` )
    | :math:`u^*` : **gU** ( :filelink:`DYNVARS.h <model/inc/DYNVARS.h>` )
    | :math:`v^*` : **gV** ( :filelink:`DYNVARS.h <model/inc/DYNVARS.h>` )
    | :math:`u^{n+1}` : **uVel** ( :filelink:`DYNVARS.h <model/inc/DYNVARS.h>` )
    | :math:`v^{n+1}` : **vVel** ( :filelink:`DYNVARS.h <model/inc/DYNVARS.h>` )


Regarding the implementation of the surface pressure solver, all
computation are done within the routine
:filelink:`SOLVE_FOR_PRESSURE <model/src/solve_for_pressure.F>` and its
dependent calls. The standard method to solve the 2D elliptic problem
:eq:`eq-solve2D` uses the conjugate gradient method
(routine :filelink:`CG2D <model/src/cg2d.F>`); the
solver matrix and conjugate gradient operator are only function of the
discretized domain and are therefore evaluated separately, before the
time iteration loop, within :filelink:`INI_CG2D <model/src/ini_cg2d.F>`. The computation of the RHS
:math:`\eta^*` is partly done in :filelink:`CALC_DIV_GHAT <model/src/calc_div_ghat.F>` and in
:filelink:`SOLVE\_FOR\_PRESSURE <model/src/solve_for_pressure.F>`.

The same method is applied for the non hydrostatic part, using a
conjugate gradient 3D solver (:filelink:`CG3D <model/src/cg3d.F>`) that is initialized in
:filelink:`INI_CG3D <model/src/ini_cg3d.F>`. The RHS terms of 2D and 3D problems are computed together
at the same point in the code.

.. toctree::
   :maxdepth:3

   crank-nicol.rst
   nonlinear-freesurf.rst


Spatial discretization of the dynamical equations
=================================================

Spatial discretization is carried out using the finite volume method.
This amounts to a grid-point method (namely second-order centered finite
difference) in the fluid interior but allows boundaries to intersect a
regular grid allowing a more accurate representation of the position of
the boundary. We treat the horizontal and vertical directions as
separable and differently.

.. toctree::
   :maxdepth:3

   finitevol-meth.rst
   c-grid.rst
   horiz-grid.rst
   vert-grid.rst
   
Continuity and horizontal pressure gradient term
=================================================


The core algorithm is based on the “C grid” discretization of the
continuity equation which can be summarized as:

.. math::
   \partial_t u + \frac{1}{\Delta x_c} \delta_i \left. \frac{ \partial \Phi}{\partial r}\right|_{s} \eta + \frac{\epsilon_{nh}}{\Delta x_c} \delta_i \Phi_{nh}' = G_u - \frac{1}{\Delta x_c} \delta_i \Phi_h'
   :label: discrete-momu

.. math::
   \partial_t v + \frac{1}{\Delta y_c} \delta_j \left. \frac{ \partial \Phi}{\partial r}\right|_{s} \eta + \frac{\epsilon_{nh}}{\Delta y_c} \delta_j \Phi_{nh}' = G_v - \frac{1}{\Delta y_c} \delta_j \Phi_h'
   :label: discrete-momv

.. math::
   \epsilon_{nh} \left( \partial_t w + \frac{1}{\Delta r_c} \delta_k \Phi_{nh}' \right) = \epsilon_{nh} G_w + \overline{b}^k - \frac{1}{\Delta r_c} \delta_k \Phi_{h}'
   :label: discrete-momw

.. math::
   \delta_i \Delta y_g \Delta r_f h_w u +
   \delta_j \Delta x_g \Delta r_f h_s v +
   \delta_k {\cal A}_c w  = {\cal A}_c \delta_k (P-E)_{r=0}
   :label: discrete-continuity

where the continuity equation has been most naturally discretized by
staggering the three components of velocity as shown in
:numref:`cgrid3d`. The grid lengths :math:`\Delta x_c` and
:math:`\Delta y_c` are the lengths between tracer points (cell centers).
The grid lengths :math:`\Delta x_g`, :math:`\Delta y_g` are the grid
lengths between cell corners. :math:`\Delta r_f` and :math:`\Delta r_c`
are the distance (in units of :math:`r`) between level interfaces
(w-level) and level centers (tracer level). The surface area presented
in the vertical is denoted :math:`{\cal
A}_c`. The factors :math:`h_w` and :math:`h_s` are non-dimensional
fractions (between 0 and 1) that represent the fraction cell depth that
is “open” for fluid flow.

The last equation, the discrete continuity equation, can be summed in
the vertical to yield the free-surface equation:

.. math::
  {\cal A}_c \partial_t \eta + \delta_i \sum_k \Delta y_g \Delta r_f h_w
   u + \delta_j \sum_k \Delta x_g \Delta r_f h_s v = {\cal
   A}_c(P-E)_{r=0}
  :label: discrete-freesurface

The source term :math:`P-E` on the rhs of continuity accounts for the
local addition of volume due to excess precipitation and run-off over
evaporation and only enters the top-level of the ocean model.

Hydrostatic balance
===================

The vertical momentum equation has the hydrostatic or quasi-hydrostatic
balance on the right hand side. This discretization guarantees that the
conversion of potential to kinetic energy as derived from the buoyancy
equation exactly matches the form derived from the pressure gradient
terms when forming the kinetic energy equation.

In the ocean, using z-coordinates, the hydrostatic balance terms are
discretized:

.. math::
   \epsilon_{nh} \partial_t w
   + g \overline{\rho'}^k + \frac{1}{\Delta z} \delta_k \Phi_h' = \ldots
   :label: discrete_hydro_ocean

In the atmosphere, using p-coordinates, hydrostatic balance is
discretized:

.. math::
   \overline{\theta'}^k + \frac{1}{\Delta \Pi} \delta_k \Phi_h' = 0
   :label: discrete_hydro_atmos

where :math:`\Delta \Pi` is the difference in Exner function between
the pressure points. The non-hydrostatic equations are not available in
the atmosphere.

The difference in approach between ocean and atmosphere occurs because
of the direct use of the ideal gas equation in forming the potential
energy conversion term :math:`\alpha \omega`. Because of the different
representation of hydrostatic balance between
ocean and atmosphere there is no elegant way to represent both systems
using an arbitrary coordinate.

The integration for hydrostatic pressure is made in the positive
:math:`r` direction (increasing k-index). For the ocean, this is from
the free-surface down and for the atmosphere this is from the ground up.

The calculations are made in the subroutine :filelink:`CALC_PHI_HYD <model/src/calc_phi_hyd.F>`. Inside
this routine, one of other of the atmospheric/oceanic form is selected
based on the string variable :varlink:`buoyancyRelation`.
