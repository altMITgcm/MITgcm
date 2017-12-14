Vector invariant momentum equations
===================================

<!– CMIREDIR:vector\_invariant\_momentum\_eqautions: –>

The finite volume method lends itself to describing the continuity and
tracer equations in curvilinear coordinate systems. However, in
curvilinear coordinates many new metric terms appear in the momentum
equations (written in Lagrangian or flux-form) making generalization far
from elegant. Fortunately, an alternative form of the equations, the
vector invariant equations are exactly that; invariant under coordinate
transformations so that they can be applied uniformly in any orthogonal
curvilinear coordinate system such as spherical coordinates, boundary
following or the conformal spherical cube system.

The non-hydrostatic vector invariant equations read:

.. math::

   \partial_t \vec{v} + ( 2\vec{\Omega} + \vec{\zeta}) \wedge \vec{v}
   - b \hat{r}
   + \vec{\nabla} B = \vec{\nabla} \cdot \vec{\bf \tau}

 which describe motions in any orthogonal curvilinear coordinate system.
Here, :math:`B` is the Bernoulli function and :math:`\vec{\zeta}=\nabla
\wedge \vec{v}` is the vorticity vector. We can take advantage of the
elegance of these equations when discretizing them and use the discrete
definitions of the grad, curl and divergence operators to satisfy
constraints. We can also consider the analogy to forming derived
equations, such as the vorticity equation, and examine how the
discretization can be adjusted to give suitable vorticity advection
among other things.

The underlying algorithm is the same as for the flux form equations. All
that has changed is the contents of the “G’s”. For the time-being, only
the hydrostatic terms have been coded but we will indicate the points
where non-hydrostatic contributions will enter:

.. math::

   \begin{aligned}
   G_u & = & G_u^{fv} + G_u^{\zeta_3 v} + G_u^{\zeta_2 w} + G_u^{\partial_x B}
   + G_u^{\partial_z \tau^x} + G_u^{h-dissip} + G_u^{v-dissip} \\
   G_v & = & G_v^{fu} + G_v^{\zeta_3 u} + G_v^{\zeta_1 w} + G_v^{\partial_y B}
   + G_v^{\partial_z \tau^y} + G_v^{h-dissip} + G_v^{v-dissip} \\
   G_w & = & G_w^{fu} + G_w^{\zeta_1 v} + G_w^{\zeta_2 u} + G_w^{\partial_z B}
   + G_w^{h-dissip} + G_w^{v-dissip}\end{aligned}

Relative vorticity
------------------

The vertical component of relative vorticity is explicitly calculated
and use in the discretization. The particular form is crucial for
numerical stability; alternative definitions break the conservation
properties of the discrete equations.

Relative vorticity is defined:

.. math::

   \zeta_3 = \frac{\Gamma}{A_\zeta}
   = \frac{1}{{\cal A}_\zeta} ( \delta_i \Delta y_c v - \delta_j \Delta x_c u )

 where :math:`{\cal A}_\zeta` is the area of the vorticity cell
presented in the vertical and :math:`\Gamma` is the circulation about
that cell.

Kinetic energy
--------------

The kinetic energy, denoted :math:`KE`, is defined:

.. math::

   KE = \frac{1}{2} ( \overline{ u^2 }^i + \overline{ v^2 }^j 
   + \epsilon_{nh} \overline{ w^2 }^k )

Coriolis terms
--------------

The potential enstrophy conserving form of the linear Coriolis terms are
written:

.. math::

   \begin{aligned}
   G_u^{fv} & = &
   \frac{1}{\Delta x_c}
   \overline{ \frac{f}{h_\zeta} }^j \overline{ \overline{ \Delta x_g h_s v }^j }^i \\
   G_v^{fu} & = & -
   \frac{1}{\Delta y_c}
   \overline{ \frac{f}{h_\zeta} }^i \overline{ \overline{ \Delta y_g h_w u }^i }^j\end{aligned}

 Here, the Coriolis parameter :math:`f` is defined at vorticity (corner)
points.

The potential enstrophy conserving form of the non-linear Coriolis terms
are written:

.. math::

   \begin{aligned}
   G_u^{\zeta_3 v} & = &
   \frac{1}{\Delta x_c}
   \overline{ \frac{\zeta_3}{h_\zeta} }^j \overline{ \overline{ \Delta x_g h_s v }^j }^i \\
   G_v^{\zeta_3 u} & = & -
   \frac{1}{\Delta y_c}
   \overline{ \frac{\zeta_3}{h_\zeta} }^i \overline{ \overline{ \Delta y_g h_w u }^i }^j\end{aligned}

The Coriolis terms can also be evaluated together and expressed in terms
of absolute vorticity :math:`f+\zeta_3`. The potential enstrophy
conserving form using the absolute vorticity is written:

.. math::

   \begin{aligned}
   G_u^{fv} + G_u^{\zeta_3 v} & = &
   \frac{1}{\Delta x_c}
   \overline{ \frac{f + \zeta_3}{h_\zeta} }^j \overline{ \overline{ \Delta x_g h_s v }^j }^i \\
   G_v^{fu} + G_v^{\zeta_3 u} & = & -
   \frac{1}{\Delta y_c}
   \overline{ \frac{f + \zeta_3}{h_\zeta} }^i \overline{ \overline{ \Delta y_g h_w u }^i }^j\end{aligned}

The distinction between using absolute vorticity or relative vorticity
is useful when constructing higher order advection schemes; monotone
advection of relative vorticity behaves differently to monotone
advection of absolute vorticity. Currently the choice of
relative/absolute vorticity, centered/upwind/high order advection is
available only through commented subroutine calls.

Shear terms
-----------

The shear terms (:math:`\zeta_2w` and :math:`\zeta_1w`) are are
discretized to guarantee that no spurious generation of kinetic energy
is possible; the horizontal gradient of Bernoulli function has to be
consistent with the vertical advection of shear:

.. math::

   \begin{aligned}
   G_u^{\zeta_2 w} & = &
   \frac{1}{ {\cal A}_w \Delta r_f h_w } \overline{
   \overline{ {\cal A}_c w }^i ( \delta_k u - \epsilon_{nh} \delta_j w )
   }^k \\
   G_v^{\zeta_1 w} & = &
   \frac{1}{ {\cal A}_s \Delta r_f h_s } \overline{
   \overline{ {\cal A}_c w }^i ( \delta_k u - \epsilon_{nh} \delta_j w )
   }^k\end{aligned}

Gradient of Bernoulli function
------------------------------

.. math::

   \begin{aligned}
   G_u^{\partial_x B} & = &
   \frac{1}{\Delta x_c} \delta_i ( \phi' + KE ) \\
   G_v^{\partial_y B} & = &
   \frac{1}{\Delta x_y} \delta_j ( \phi' + KE )
   %G_w^{\partial_z B} & = &
   %\frac{1}{\Delta r_c} h_c \delta_k ( \phi' + KE )\end{aligned}

Horizontal divergence
---------------------

The horizontal divergence, a complimentary quantity to relative
vorticity, is used in parameterizing the Reynolds stresses and is
discretized:

.. math::

   D = \frac{1}{{\cal A}_c h_c} (
     \delta_i \Delta y_g h_w u
   + \delta_j \Delta x_g h_s v )

Horizontal dissipation
----------------------

The following discretization of horizontal dissipation conserves
potential vorticity (thickness weighted relative vorticity) and
divergence and dissipates energy, enstrophy and divergence squared:

.. math::

   \begin{aligned}
   G_u^{h-dissip} & = &
     \frac{1}{\Delta x_c} \delta_i ( A_D D - A_{D4} D^*)
   - \frac{1}{\Delta y_u h_w} \delta_j h_\zeta ( A_\zeta \zeta - A_{\zeta4} \zeta^* )
   \\
   G_v^{h-dissip} & = &
     \frac{1}{\Delta x_v h_s} \delta_i h_\zeta ( A_\zeta \zeta - A_\zeta \zeta^* )
   + \frac{1}{\Delta y_c} \delta_j ( A_D D - A_{D4} D^* )\end{aligned}

 where

.. math::

   \begin{aligned}
   D^* & = & \frac{1}{{\cal A}_c h_c} (
     \delta_i \Delta y_g h_w \nabla^2 u
   + \delta_j \Delta x_g h_s \nabla^2 v ) \\
   \zeta^* & = & \frac{1}{{\cal A}_\zeta} (
     \delta_i \Delta y_c \nabla^2 v
   - \delta_j \Delta x_c \nabla^2 u )\end{aligned}

Vertical dissipation
--------------------

Currently, this is exactly the same code as the flux form equations.

.. math::

   \begin{aligned}
   G_u^{v-diss} & = &
   \frac{1}{\Delta r_f h_w} \delta_k \tau_{13} \\
   G_v^{v-diss} & = &
   \frac{1}{\Delta r_f h_s} \delta_k \tau_{23}\end{aligned}

 represents the general discrete form of the vertical dissipation terms.

In the interior the vertical stresses are discretized:

.. math::

   \begin{aligned}
   \tau_{13} & = & A_v \frac{1}{\Delta r_c} \delta_k u \\
   \tau_{23} & = & A_v \frac{1}{\Delta r_c} \delta_k v\end{aligned}
