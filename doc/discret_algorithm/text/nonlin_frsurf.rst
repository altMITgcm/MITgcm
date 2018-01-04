Non-linear free-surface
-----------------------

Recently, new options have been added to the model that concern the free
surface formulation.

pressure/geo-potential and free surface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| For the atmosphere, since
  :math:`\phi = \phi_{topo} - \int^p_{p_s} \alpha dp`, subtracting the
  reference state defined in section [sec:hpe-p-geo-potential-split] :
| 

  .. math::

     \phi_o = \phi_{topo} - \int^p_{p_o} \alpha_o dp
     \hspace{5mm}\mathrm{with}\hspace{3mm} \phi_o(p_o)=\phi_{topo}

   we get:

  .. math:: \phi' = \phi - \phi_o = \int^{p_s}_p \alpha dp - \int^{p_o}_p \alpha_o dp

   For the ocean, the reference state is simpler since :math:`\rho_c`
  does not dependent on :math:`z` (:math:`b_o=g`) and the surface
  reference position is uniformly :math:`z=0` (:math:`R_o=0`), and the
  same subtraction leads to a similar relation. For both fluid, using
  the isomorphic notations, we can write:

  .. math:: \phi' = \int^{r_{surf}}_r b~ dr - \int^{R_o}_r b_o dr

   and re-write as:

  .. math::

     \phi' = \int^{r_{surf}}_{R_o} b~ dr + \int^{R_o}_r (b - b_o) dr
     \label{eq:split-phi-Ro}

   or:

  .. math::

     \phi' = \int^{r_{surf}}_{R_o} b_o dr + \int^{r_{surf}}_r (b - b_o) dr
     \label{eq:split-phi-bo}

In section [sec:finding\_the\_pressure\_field], following
eq.[eq:split-phi-Ro], the pressure/geo-potential :math:`\phi'` has been
separated into surface (:math:`\phi_s`), and hydrostatic anomaly
(:math:`\phi'_{hyd}`). In this section, the split between :math:`\phi_s`
and :math:`\phi'_{hyd}` is made according to equation [eq:split-phi-bo].
This slightly different definition reflects the actual implementation in
the code and is valid for both linear and non-linear free-surface
formulation, in both r-coordinate and r\*-coordinate.

| Because the linear free-surface approximation ignore the tracer
  content of the fluid parcel between :math:`R_o` and
  :math:`r_{surf}=R_o+\eta`, for consistency reasons, this part is also
  neglected in :math:`\phi'_{hyd}` :

  .. math:: \phi'_{hyd} = \int^{r_{surf}}_r (b - b_o) dr \simeq \int^{R_o}_r (b - b_o) dr

   Note that in this case, the two definitions of :math:`\phi_s` and
  :math:`\phi'_{hyd}` from equation [eq:split-phi-Ro] and
  [eq:split-phi-bo] converge toward the same (approximated) expressions:
  :math:`\phi_s = \int^{r_{surf}}_{R_o} b_o dr` and
  :math:`\phi'_{hyd}=\int^{R_o}_r b' dr`.
| On the contrary, the unapproximated formulation (“non-linear
  free-surface”, see the next section) retains the full expression:
  :math:`\phi'_{hyd} = \int^{r_{surf}}_r (b - b_o) dr ` . This is
  obtained by selecting **nonlinFreeSurf**\ =4 in parameter file *data*.
| Regarding the surface potential:

  .. math::

     \phi_s = \int_{R_o}^{R_o+\eta} b_o dr = b_s \eta
     \hspace{5mm}\mathrm{with}\hspace{5mm}
     b_s = \frac{1}{\eta} \int_{R_o}^{R_o+\eta} b_o dr

   :math:`b_s \simeq b_o(R_o)` is an excellent approximation (better
  than the usual numerical truncation, since generally :math:`|\eta|` is
  smaller than the vertical grid increment).

For the ocean, :math:`\phi_s = g \eta` and :math:`b_s = g` is uniform.
For the atmosphere, however, because of topographic effects, the
reference surface pressure :math:`R_o=p_o` has large spatial variations
that are responsible for significant :math:`b_s` variations (from 0.8 to
1.2 :math:`[m^3/kg]`). For this reason, when **uniformLin\_PhiSurf**
*=.FALSE.* (parameter file *data*, namelist *PARAM01*) a non-uniform
linear coefficient :math:`b_s` is used and computed (*S/R
INI\_LINEAR\_PHISURF*) according to the reference surface pressure
:math:`p_o`:
:math:`b_s = b_o(R_o) = c_p \kappa (p_o / P^o_{SL})^{(\kappa - 1)} \theta_{ref}(p_o)`.
with :math:`P^o_{SL}` the mean sea-level pressure.

Free surface effect on column total thickness (Non-linear free-surface)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The total thickness of the fluid column is :math:`r_{surf} - R_{fixed} =
\eta + R_o - R_{fixed}`. In most applications, the free surface
displacements are small compared to the total thickness
:math:`\eta \ll H_o = R_o - R_{fixed}`. In the previous sections and in
older version of the model, the linearized free-surface approximation
was made, assuming :math:`r_{surf} - R_{fixed} \simeq H_o` when
computing horizontal transports, either in the continuity equation or in
tracer and momentum advection terms. This approximation is dropped when
using the non-linear free-surface formulation and the total thickness,
including the time varying part :math:`\eta`, is considered when
computing horizontal transports. Implications for the barotropic part
are presented hereafter. In section [sec:freesurf-tracer-advection]
consequences for tracer conservation is briefly discussed (more details
can be found in :raw-latex:`\cite{campin:02}`) ; the general
time-stepping is presented in section [sec:nonlin-freesurf-timestepping]
with some limitations regarding the vertical resolution in section
[sec:nonlin-freesurf-dz\_surf].

In the non-linear formulation, the continuous form of the model
equations remains unchanged, except for the 2D continuity equation
([eq:discrete-time-backward-free-surface]) which is now integrated from
:math:`R_{fixed}(x,y)` up to :math:`r_{surf}=R_o+\eta` :

.. math::

   \epsilon_{fs} \partial_t \eta =
   \left. \dot{r} \right|_{r=r_{surf}} + \epsilon_{fw} (P-E) =
   - {\bf \nabla}_h \cdot \int_{R_{fixed}}^{R_o+\eta} \vec{\bf v} dr
   + \epsilon_{fw} (P-E)

Since :math:`\eta` has a direct effect on the horizontal velocity
(through :math:`\nabla_h \Phi_{surf}`), this adds a non-linear term to
the free surface equation. Several options for the time discretization
of this non-linear part can be considered, as detailed below.

If the column thickness is evaluated at time step :math:`n`, and with
implicit treatment of the surface potential gradient, equations
([eq-solve2D] and [eq-solve2D\_rhs]) becomes:

.. math::

   \begin{aligned}
   \epsilon_{fs} {\eta}^{n+1} -
   {\bf \nabla}_h \cdot \Delta t^2 (\eta^{n}+R_o-R_{fixed})
   {\bf \nabla}_h b_s {\eta}^{n+1}
   = {\eta}^*\end{aligned}

 where

.. math::

   \begin{aligned}
   {\eta}^* = \epsilon_{fs} \: {\eta}^{n} -
   \Delta t {\bf \nabla}_h \cdot \int_{R_{fixed}}^{R_o+\eta^n} \vec{\bf v}^* dr
   \: + \: \epsilon_{fw} \Delta_t (P-E)^{n}\end{aligned}

 This method requires us to update the solver matrix at each time step.

Alternatively, the non-linear contribution can be evaluated fully
explicitly:

.. math::

   \begin{aligned}
   \epsilon_{fs} {\eta}^{n+1} -
   {\bf \nabla}_h \cdot \Delta t^2 (R_o-R_{fixed})
   {\bf \nabla}_h b_s {\eta}^{n+1}
   = {\eta}^*
   +{\bf \nabla}_h \cdot \Delta t^2 (\eta^{n})
   {\bf \nabla}_h b_s {\eta}^{n}\end{aligned}

 This formulation allows one to keep the initial solver matrix unchanged
though throughout the integration, since the non-linear free surface
only affects the RHS.

Finally, another option is a “linearized” formulation where the total
column thickness appears only in the integral term of the RHS
([eq-solve2D\_rhs]) but not directly in the equation ([eq-solve2D]).

Those different options (see Table [tab:nonLinFreeSurf\_flags]) have
been tested and show little differences. However, we recommend the use
of the most precise method (the 1rst one) since the computation cost
involved in the solver matrix update is negligible.

+------------------+---------+----------------------------------------------------------------------------------------+
| parameter        | value   | description                                                                            |
+==================+=========+========================================================================================+
|                  | -1      | linear free-surface, restart from a pickup file                                        |
+------------------+---------+----------------------------------------------------------------------------------------+
|                  |         | produced with #undef EXACT\_CONSERV code                                               |
+------------------+---------+----------------------------------------------------------------------------------------+
|                  | 0       | Linear free-surface                                                                    |
+------------------+---------+----------------------------------------------------------------------------------------+
| nonlinFreeSurf   | 4       | Non-linear free-surface                                                                |
+------------------+---------+----------------------------------------------------------------------------------------+
|                  | 3       | same as 4 but neglecting :math:`\int_{R_o}^{R_o+\eta} b' dr ` in :math:`\Phi'_{hyd}`   |
+------------------+---------+----------------------------------------------------------------------------------------+
|                  | 2       | same as 3 but do not update cg2d solver matrix                                         |
+------------------+---------+----------------------------------------------------------------------------------------+
|                  | 1       | same as 2 but treat momentum as in Linear FS                                           |
+------------------+---------+----------------------------------------------------------------------------------------+
|                  | 0       | do not use :math:`r*` vertical coordinate (= default)                                  |
+------------------+---------+----------------------------------------------------------------------------------------+
| select\_rStar    | 2       | use :math:`r^*` vertical coordinate                                                    |
+------------------+---------+----------------------------------------------------------------------------------------+
|                  | 1       | same as 2 but without the contribution of the                                          |
+------------------+---------+----------------------------------------------------------------------------------------+
|                  |         | slope of the coordinate in :math:`\nabla \Phi`                                         |
+------------------+---------+----------------------------------------------------------------------------------------+

Table: Non-linear free-surface flags

Tracer conservation with non-linear free-surface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To ensure global tracer conservation (i.e., the total amount) as well as
local conservation, the change in the surface level thickness must be
consistent with the way the continuity equation is integrated, both in
the barotropic part (to find :math:`\eta`) and baroclinic part (to find
:math:`w = \dot{r}`).

To illustrate this, consider the shallow water model, with a source of
fresh water (P):

.. math:: \partial_t h + \nabla \cdot h \vec{\bf v} = P

 where :math:`h` is the total thickness of the water column. To conserve
the tracer :math:`\theta` we have to discretize:

.. math::

   \partial_t (h \theta) + \nabla \cdot ( h \theta \vec{\bf v})
     = P \theta_{\mathrm{rain}}

 Using the implicit (non-linear) free surface described above (section
[sec:pressure-method-linear-backward]) we have:

.. math::

   \begin{aligned}
   h^{n+1} = h^{n} - \Delta t \nabla \cdot (h^n \, \vec{\bf v}^{n+1} ) + \Delta t P \\\end{aligned}

 The discretized form of the tracer equation must adopt the same “form”
in the computation of tracer fluxes, that is, the same value of
:math:`h`, as used in the continuity equation:

.. math::

   \begin{aligned}
   h^{n+1} \, \theta^{n+1} = h^n \, \theta^n
           - \Delta t \nabla \cdot (h^n \, \theta^n \, \vec{\bf v}^{n+1})
           + \Delta t P \theta_{rain}\end{aligned}

The use of a 3 time-levels time-stepping scheme such as the
Adams-Bashforth make the conservation sightly tricky. The current
implementation with the Adams-Bashforth time-stepping provides an exact
local conservation and prevents any drift in the global tracer content
(:raw-latex:`\cite{campin:02}`). Compared to the linear free-surface
method, an additional step is required: the variation of the water
column thickness (from :math:`h^n` to :math:`h^{n+1}`) is not
incorporated directly into the tracer equation. Instead, the model uses
the :math:`G_\theta` terms (first step) as in the linear free surface
formulation (with the “*surface correction*” turned “on”, see tracer
section):

.. math::

   G_\theta^n = \left(- \nabla \cdot (h^n \, \theta^n \, \vec{\bf v}^{n+1})
            - \dot{r}_{surf}^{n+1} \theta^n \right) / h^n

 Then, in a second step, the thickness variation (expansion/reduction)
is taken into account:

.. math::

   \theta^{n+1} = \theta^n + \Delta t \frac{h^n}{h^{n+1}}
      \left( G_\theta^{(n+1/2)} + P (\theta_{\mathrm{rain}} - \theta^n )/h^n \right)
   %\theta^{n+1} = \theta^n + \frac{\Delta t}{h^{n+1}}
   %   \left( h^n G_\theta^{(n+1/2)} + P (\theta_{\mathrm{rain}} - \theta^n ) \right)

 Note that with a simple forward time step (no Adams-Bashforth), these
two formulations are equivalent, since :math:`
(h^{n+1} - h^{n})/ \Delta t =
P - \nabla \cdot (h^n \, \vec{\bf v}^{n+1} ) = P + \dot{r}_{surf}^{n+1}
`

Time stepping implementation of the non-linear free-surface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The grid cell thickness was hold constant with the linear free-surface ;
with the non-linear free-surface, it is now varying in time, at least at
the surface level. This implies some modifications of the general
algorithm described earlier in sections [sec:adams-bashforth-sync] and
[sec:adams-bashforth-staggered].

| A simplified version of the staggered in time, non-linear free-surface
  algorithm is detailed hereafter, and can be compared to the equivalent
  linear free-surface case (eq. [eq:Gv-n-staggered] to
  [eq:t-n+1-staggered]) and can also be easily transposed to the
  synchronous time-stepping case. Among the simplifications, salinity
  equation, implicit operator and detailed elliptic equation are
  omitted. Surface forcing is explicitly written as fluxes of
  temperature, fresh water and momentum,
  :math:`Q^{n+1/2}, P^{n+1/2}, F_{\bf v}^n` respectively. :math:`h^n`
  and :math:`dh^n` are the column and grid box thickness in
  r-coordinate.

  .. math::

     \begin{aligned}
     \phi^{n}_{hyd} & = & \int b(\theta^{n},S^{n},r) dr
     \label{eq:phi-hyd-nlfs} \\
     \vec{\bf G}_{\vec{\bf v}}^{n-1/2}\hspace{-2mm} & = &
     \vec{\bf G}_{\vec{\bf v}} (dh^{n-1},\vec{\bf v}^{n-1/2})
     \hspace{+2mm};\hspace{+2mm}
     \vec{\bf G}_{\vec{\bf v}}^{(n)} =
        \frac{3}{2} \vec{\bf G}_{\vec{\bf v}}^{n-1/2}
     -  \frac{1}{2} \vec{\bf G}_{\vec{\bf v}}^{n-3/2}
     \label{eq:Gv-n-nlfs} \\
     %\vec{\bf G}_{\vec{\bf v}}^{(n)} & = &
     %   \frac{3}{2} \vec{\bf G}_{\vec{\bf v}}^{n-1/2}
     %-  \frac{1}{2} \vec{\bf G}_{\vec{\bf v}}^{n-3/2}
     %\label{eq:Gv-n+5-nlfs} \\
     %\vec{\bf v}^{*} & = & \vec{\bf v}^{n-1/2} + \frac{\Delta t}{dh^{n}} \left(
     %dh^{n-1}\vec{\bf G}_{\vec{\bf v}}^{(n)} + F_{\vec{\bf v}}^{n} \right)
     \vec{\bf v}^{*} & = & \vec{\bf v}^{n-1/2} + \Delta t \frac{dh^{n-1}}{dh^{n}} \left(
     \vec{\bf G}_{\vec{\bf v}}^{(n)} + F_{\vec{\bf v}}^{n}/dh^{n-1} \right)
     - \Delta t \nabla \phi_{hyd}^{n}
     \label{eq:vstar-nlfs}\end{aligned}

   :math:`\longrightarrow`  *update
  model geometry : *\ :math:`{\bf hFac}(dh^n)`
| 

  .. math::

     \begin{aligned}
     \eta^{n+1/2} \hspace{-2mm} & = &
     \eta^{n-1/2} + \Delta t P^{n+1/2} - \Delta t
       \nabla \cdot \int \vec{\bf v}^{n+1/2} dh^{n} \nonumber \\
                  & = & \eta^{n-1/2} + \Delta t P^{n+1/2} - \Delta t
       \nabla \cdot \int \!\!\! \left( \vec{\bf v}^* - g \Delta t \nabla \eta^{n+1/2} \right) dh^{n}
     \label{eq:nstar-nlfs} \\
     \vec{\bf v}^{n+1/2}\hspace{-2mm} & = &
     \vec{\bf v}^{*} - g \Delta t \nabla \eta^{n+1/2}
     \label{eq:v-n+1-nlfs} \\
     h^{n+1} & = & h^{n} + \Delta t P^{n+1/2} - \Delta t
       \nabla \cdot \int \vec{\bf v}^{n+1/2} dh^{n}
     \label{eq:h-n+1-nlfs} \\
     G_{\theta}^{n} & = & G_{\theta} ( dh^{n}, u^{n+1/2}, \theta^{n} )
     \hspace{+2mm};\hspace{+2mm}
     G_{\theta}^{(n+1/2)} = \frac{3}{2} G_{\theta}^{n} - \frac{1}{2} G_{\theta}^{n-1}
     \label{eq:Gt-n-nlfs} \\
     %\theta^{n+1} & = &\theta^{n} + \frac{\Delta t}{dh^{n+1}} \left( dh^n
     %G_{\theta}^{(n+1/2)} + Q^{n+1/2} + P^{n+1/2} (\theta_{\mathrm{rain}}-\theta^n) \right)
     \theta^{n+1} & = &\theta^{n} + \Delta t \frac{dh^n}{dh^{n+1}} \left(
     G_{\theta}^{(n+1/2)}
     +( P^{n+1/2} (\theta_{\mathrm{rain}}-\theta^n) + Q^{n+1/2})/dh^n \right)
     \nonumber \\
     & & \label{eq:t-n+1-nlfs}\end{aligned}

   Two steps have been added to linear free-surface algorithm (eq.
  [eq:Gv-n-staggered] to [eq:t-n+1-staggered]): Firstly, the model
  “geometry” (here the **hFacC,W,S**) is updated just before entering
  *SOLVE\_FOR\_PRESSURE*, using the current :math:`dh^{n}` field.
  Secondly, the vertically integrated continuity equation
  (eq.[eq:h-n+1-nlfs]) has been added (**exactConserv**\ *=TRUE*, in
  parameter file *data*, namelist *PARM01*) just before computing the
  vertical velocity, in subroutine *INTEGR\_CONTINUITY*. Although this
  equation might appear redundant with eq.[eq:nstar-nlfs], the
  integrated column thickness :math:`h^{n+1}` will be different from
  :math:`\eta^{n+1/2} + H`  in the following cases:

-  when Crank-Nicolson time-stepping is used (see section
   [sec:freesurf-CrankNick]).

-  when filters are applied to the flow field, after ([eq:v-n+1-nlfs])
   and alter the divergence of the flow.

-  when the solver does not iterate until convergence ; for example,
   because a too large residual target was set (**cg2dTargetResidual**,
   parameter file *data*, namelist *PARM02*).

In this staggered time-stepping algorithm, the momentum tendencies are
computed using :math:`dh^{n-1}` geometry factors. (eq.[eq:Gv-n-nlfs])
and then rescaled in subroutine *TIMESTEP*, (eq.[eq:vstar-nlfs]),
similarly to tracer tendencies (see section
[sec:freesurf-tracer-advection]). The tracers are stepped forward later,
using the recently updated flow field :math:`{\bf v}^{n+1/2}` and the
corresponding model geometry :math:`dh^{n}` to compute the tendencies
(eq.[eq:Gt-n-nlfs]); Then the tendencies are rescaled by
:math:`dh^n/dh^{n+1}` to derive the new tracers values
:math:`(\theta,S)^{n+1}` (eq.[eq:t-n+1-nlfs], in subroutine *CALC\_GT,
CALC\_GS*).

Note that the fresh-water input is added in a consistent way in the
continuity equation and in the tracer equation, taking into account the
fresh-water temperature :math:`\theta_{\mathrm{rain}}`.

| Regarding the restart procedure, two 2.D fields :math:`h^{n-1}` and
  :math:`(h^n-h^{n-1})/\Delta t` in addition to the standard state
  variables and tendencies (:math:`\eta^{n-1/2}`,
  :math:`{\bf v}^{n-1/2}`, :math:`\theta^n`, :math:`S^n`,
  :math:`{\bf G}_{\bf v}^{n-3/2}`, :math:`G_{\theta,S}^{n-1}`) are
  stored in a “*pickup*” file. The model restarts reading this
  “*pickup*” file, then update the model geometry according to
  :math:`h^{n-1}`, and compute :math:`h^n` and the vertical velocity
  before starting the main calling sequence (eq.[eq:phi-hyd-nlfs] to
  [eq:t-n+1-nlfs], *S/R FORWARD\_STEP*).

Non-linear free-surface and vertical resolution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the amplitude of the free-surface variations becomes as large as
the vertical resolution near the surface, the surface layer thickness
can decrease to nearly zero or can even vanish completely. This later
possibility has not been implemented, and a minimum relative thickness
is imposed (**hFacInf**, parameter file *data*, namelist *PARM01*) to
prevent numerical instabilities caused by very thin surface level.

A better alternative to the vanishing level problem has been found and
implemented recently, and rely on a different vertical coordinate
:math:`r^*` : The time variation ot the total column thickness becomes
part of the r\* coordinate motion, as in a :math:`\sigma_{z},\sigma_{p}`
model, but the fixed part related to topography is treated as in a
height or pressure coordinate model. A complete description is given in
:raw-latex:`\cite{adcroft:04a}`.

The time-stepping implementation of the :math:`r^*` coordinate is
identical to the non-linear free-surface in :math:`r` coordinate, and
differences appear only in the spacial discretization.
