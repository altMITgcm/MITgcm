Continuous equations in ‘r’ coordinates
=======================================
To render atmosphere and ocean models from one dynamical core we exploit
‘isomorphisms’ between equation sets that govern the evolution of the
respective fluids - see :numref:`isomorphic-equations`. One system of
hydrodynamical equations is written down and encoded. The model
variables have different interpretations depending on whether the
atmosphere or ocean is being studied. Thus, for example, the vertical
coordinate ‘:math:`r`’ is interpreted as pressure, :math:`p`, if we are
modeling the atmosphere (right hand side of :numref:`isomorphic-equations`) and height, :math:`z`, if we are modeling
the ocean (left hand side of :numref:`isomorphic-equations`).


  .. figure:: figs/zandpcoord.*
    :width: 90%
    :align: center
    :alt: isomorphic-equations
    :name: isomorphic-equations

    Isomorphic equation sets used for atmosphere (right) and ocean (left).


The state of the fluid at any time is characterized by the distribution
of velocity :math:`\vec{\mathbf{v}}`, active tracers :math:`\theta` and
:math:`S`, a ‘geopotential’ :math:`\phi` and density
:math:`\rho =\rho (\theta ,S,p)` which may depend on :math:`\theta`,
:math:`S`, and :math:`p`. The equations that govern the evolution of
these fields, obtained by applying the laws of classical mechanics and
thermodynamics to a Boussinesq, Navier-Stokes fluid are, written in
terms of a generic vertical coordinate, :math:`r`, so that the
appropriate kinematic boundary conditions can be applied isomorphically
see :numref:`zandp-vert-coord`.


  .. figure:: figs/vertcoord.*
    :width: 70%
    :align: center
    :alt: zandp-vert-coord
    :name: zandp-vert-coord

    Vertical coordinates and kinematic boundary conditions for atmosphere (top) and ocean (bottom).

.. math::
   \frac{D\vec{\mathbf{v}_{h}}}{Dt}+\left( 2\vec{\Omega}\times \vec{\mathbf{v}}
   \right) _{h}+\mathbf{\nabla }_{h}\phi =\mathcal{F}_{\vec{\mathbf{v}_{h}}}\text{  horizontal momentum}
   :label: horiz-mtm

.. math::
   \frac{D\dot{r}}{Dt}+\widehat{k}\cdot \left( 2\vec{\Omega}\times \vec{\mathbf{
   v}}\right) +\frac{\partial \phi }{\partial r}+b=\mathcal{F}_{\dot{r}}\text{  vertical momentum}
   :label: vert-mtm

.. math::
   \mathbf{\nabla }_{h}\cdot \vec{\mathbf{v}}_{h}+\frac{\partial \dot{r}}{
   \partial r}=0\text{  continuity}
   :label: continuity

.. math:: 
   b=b(\theta ,S,r)\text{  equation of state} 
   :label: eos
 
.. math::
   \frac{D\theta }{Dt}=\mathcal{Q}_{\theta }\text{  potential temperature}
   :label: pot-temp

.. math::
   \frac{DS}{Dt}=\mathcal{Q}_{S}\text{  humidity/salinity}
   :label: humidity-salt

Here:

.. math:: r\text{ is the vertical coordinate}

.. math::

   \frac{D}{Dt}=\frac{\partial }{\partial t}+\vec{\mathbf{v}}\cdot \nabla \text{ is the total derivative}

.. math::

   \mathbf{\nabla }=\mathbf{\nabla }_{h}+\widehat{k}\frac{\partial }{\partial r}
   \text{  is the ‘grad’ operator}

with :math:`\mathbf{\nabla }_{h}` operating in the horizontal and
:math:`\widehat{k}
\frac{\partial }{\partial r}` operating in the vertical, where
:math:`\widehat{k}` is a unit vector in the vertical

.. math:: t\text{ is time}

.. math::

   \vec{\mathbf{v}}=(u,v,\dot{r})=(\vec{\mathbf{v}}_{h},\dot{r})\text{ is the velocity}

.. math:: \phi \text{ is the ‘pressure’/‘geopotential’}

.. math:: \vec{\Omega}\text{ is the Earth's rotation}

.. math:: b\text{ is the ‘buoyancy’}

.. math:: \theta \text{ is potential temperature}

.. math:: S\text{ is specific humidity in the atmosphere; salinity in the ocean}

.. math::

   \mathcal{F}_{\vec{\mathbf{v}}}\text{ are forcing and dissipation of }\vec{
   \mathbf{v}}

.. math:: \mathcal{Q}_{\theta }\mathcal{\ }\text{ are forcing and dissipation of }\theta

.. math:: \mathcal{Q}_{S}\mathcal{\ }\text{are forcing and dissipation of }S

The :math:`\mathcal{F}^{\prime }s` and :math:`\mathcal{Q}^{\prime }s`
are provided by ‘physics’ and forcing packages for atmosphere and ocean.
These are described in later chapters.


.. toctree::
   :maxdepth: 3

   kinematic_bound.rst
   atmosphere.rst
   ocean.rst
   hydrostatic.rst
   soln_strategy.rst
   finding_pressure.rst
   forcing_dissip.rst
   vector_invar.rst
   adjoint.rst