Spatial discretization of the dynamical equations
=================================================

<!– CMIREDIR:spatial\_discretization\_of\_dyn\_eq: –>

Spatial discretization is carried out using the finite volume method.
This amounts to a grid-point method (namely second-order centered finite
difference) in the fluid interior but allows boundaries to intersect a
regular grid allowing a more accurate representation of the position of
the boundary. We treat the horizontal and vertical directions as
separable and differently.

The finite volume method: finite volumes versus finite difference
-----------------------------------------------------------------

<!– CMIREDIR:finite\_volume: –>

The finite volume method is used to discretize the equations in space.
The expression “finite volume” actually has two meanings; one is the
method of embedded or intersecting boundaries (shaved or lopped cells in
our terminology) and the other is non-linear interpolation methods that
can deal with non-smooth solutions such as shocks (i.e. flux limiters
for advection). Both make use of the integral form of the conservation
laws to which the *weak solution* is a solution on each finite volume of
(sub-domain). The weak solution can be constructed out of piece-wise
constant elements or be differentiable. The differentiable equations can
not be satisfied by piece-wise constant functions.

As an example, the 1-D constant coefficient advection-diffusion
equation:

.. math:: \partial_t \theta + \partial_x ( u \theta - \kappa \partial_x \theta ) = 0

 can be discretized by integrating over finite sub-domains, i.e. the
lengths :math:`\Delta x_i`:

.. math:: \Delta x \partial_t \theta + \delta_i ( F ) = 0

 is exact if :math:`\theta(x)` is piece-wise constant over the interval
:math:`\Delta x_i` or more generally if :math:`\theta_i` is defined as
the average over the interval :math:`\Delta x_i`.

The flux, :math:`F_{i-1/2}`, must be approximated:

.. math:: F = u \overline{\theta} - \frac{\kappa}{\Delta x_c} \partial_i \theta

 and this is where truncation errors can enter the solution. The method
for obtaining :math:`\overline{\theta}` is unspecified and a wide range
of possibilities exist including centered and upwind interpolation,
polynomial fits based on the the volume average definitions of
quantities and non-linear interpolation such as flux-limiters.

Choosing simple centered second-order interpolation and differencing
recovers the same ODE’s resulting from finite differencing for the
interior of a fluid. Differences arise at boundaries where a boundary is
not positioned on a regular or smoothly varying grid. This method is
used to represent the topography using lopped cell, see
:raw-latex:`\cite{adcroft:97}`. Subtle difference also appear in more
than one dimension away from boundaries. This happens because the each
direction is discretized independently in the finite difference method
while the integrating over finite volume implicitly treats all
directions simultaneously. Illustration of this is given in
:raw-latex:`\cite{ac:02}`.

C grid staggering of variables
------------------------------

The basic algorithm employed for stepping forward the momentum equations
is based on retaining non-divergence of the flow at all times. This is
most naturally done if the components of flow are staggered in space in
the form of an Arakawa C grid :raw-latex:`\cite{arakawa:77}`.

Fig. [fig:cgrid3d] shows the components of flow
(:math:`u`,\ :math:`v`,\ :math:`w`) staggered in space such that the
zonal component falls on the interface between continuity cells in the
zonal direction. Similarly for the meridional and vertical directions.
The continuity cell is synonymous with tracer cells (they are one and
the same).

Grid initialization and data
----------------------------

Initialization of grid data is controlled by subroutine *INI\_GRID*
which in calls *INI\_VERTICAL\_GRID* to initialize the vertical grid,
and then either of *INI\_CARTESIAN\_GRID*, *INI\_SPHERICAL\_POLAR\_GRID*
or *INI\_CURVILINEAR\_GRID* to initialize the horizontal grid for
cartesian, spherical-polar or curvilinear coordinates respectively.

The reciprocals of all grid quantities are pre-calculated and this is
done in subroutine *INI\_MASKS\_ETC* which is called later by subroutine
*INITIALIZE\_FIXED*.

All grid descriptors are global arrays and stored in common blocks in
*GRID.h* and a generally declared as *\_RS*.

Horizontal grid
---------------

+----+----+
+----+----+
+----+----+

The model domain is decomposed into tiles and within each tile a
quasi-regular grid is used. A tile is the basic unit of domain
decomposition for parallelization but may be used whether parallelized
or not; see section [sec:domain\_decomposition] for more details.
Although the tiles may be patched together in an unstructured manner
(i.e. irregular or non-tessilating pattern), the interior of tiles is a
structured grid of quadrilateral cells. The horizontal coordinate system
is orthogonal curvilinear meaning we can not necessarily treat the two
horizontal directions as separable. Instead, each cell in the horizontal
grid is described by the length of it’s sides and it’s area.

The grid information is quite general and describes any of the available
coordinates systems, cartesian, spherical-polar or curvilinear. All that
is necessary to distinguish between the coordinate systems is to
initialize the grid data (descriptors) appropriately.

In the following, we refer to the orientation of quantities on the
computational grid using geographic terminology such as points of the
compass. This is purely for convenience but should not be confused with
the actual geographic orientation of model quantities.

Fig. [fig:hgrid]a shows the tracer cell (synonymous with the continuity
cell). The length of the southern edge, :math:`\Delta x_g`, western
edge, :math:`\Delta y_g` and surface area, :math:`A_c`, presented in the
vertical are stored in arrays **DXg**, **DYg** and **rAc**. The “g”
suffix indicates that the lengths are along the defining grid
boundaries. The “c” suffix associates the quantity with the cell
centers. The quantities are staggered in space and the indexing is such
that **DXg(i,j)** is positioned to the south of **rAc(i,j)** and
**DYg(i,j)** positioned to the west.

Fig. [fig:hgrid]b shows the vorticity cell. The length of the southern
edge, :math:`\Delta x_c`, western edge, :math:`\Delta y_c` and surface
area, :math:`A_\zeta`, presented in the vertical are stored in arrays
**DXc**, **DYc** and **rAz**. The “z” suffix indicates that the lengths
are measured between the cell centers and the “:math:`\zeta`” suffix
associates points with the vorticity points. The quantities are
staggered in space and the indexing is such that **DXc(i,j)** is
positioned to the north of **rAz(i,j)** and **DYc(i,j)** positioned to
the east.

Fig. [fig:hgrid]c shows the “u” or western (w) cell. The length of the
southern edge, :math:`\Delta x_v`, eastern edge, :math:`\Delta y_f` and
surface area, :math:`A_w`, presented in the vertical are stored in
arrays **DXv**, **DYf** and **rAw**. The “v” suffix indicates that the
length is measured between the v-points, the “f” suffix indicates that
the length is measured between the (tracer) cell faces and the “w”
suffix associates points with the u-points (w stands for west). The
quantities are staggered in space and the indexing is such that
**DXv(i,j)** is positioned to the south of **rAw(i,j)** and **DYf(i,j)**
positioned to the east.

Fig. [fig:hgrid]d shows the “v” or southern (s) cell. The length of the
northern edge, :math:`\Delta x_f`, western edge, :math:`\Delta y_u` and
surface area, :math:`A_s`, presented in the vertical are stored in
arrays **DXf**, **DYu** and **rAs**. The “u” suffix indicates that the
length is measured between the u-points, the “f” suffix indicates that
the length is measured between the (tracer) cell faces and the “s”
suffix associates points with the v-points (s stands for south). The
quantities are staggered in space and the indexing is such that
**DXf(i,j)** is positioned to the north of **rAs(i,j)** and **DYu(i,j)**
positioned to the west.

Reciprocals of horizontal grid descriptors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lengths and areas appear in the denominator of expressions as much as in
the numerator. For efficiency and portability, we pre-calculate the
reciprocal of the horizontal grid quantities so that in-line divisions
can be avoided.

For each grid descriptor (array) there is a reciprocal named using the
prefix **RECIP\_**. This doubles the amount of storage in *GRID.h* but
they are all only 2-D descriptors.

Cartesian coordinates
~~~~~~~~~~~~~~~~~~~~~

Cartesian coordinates are selected when the logical flag
**usingCartesianGrid** in namelist *PARM04* is set to true. The grid
spacing can be set to uniform via scalars **dXspacing** and
**dYspacing** in namelist *PARM04* or to variable resolution by the
vectors **DELX** and **DELY**. Units are normally meters.
Non-dimensional coordinates can be used by interpreting the
gravitational constant as the Rayleigh number.

Spherical-polar coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spherical coordinates are selected when the logical flag
**usingSphericalPolarGrid** in namelist *PARM04* is set to true. The
grid spacing can be set to uniform via scalars **dXspacing** and
**dYspacing** in namelist *PARM04* or to variable resolution by the
vectors **DELX** and **DELY**. Units of these namelist variables are
alway degrees. The horizontal grid descriptors are calculated from these
namelist variables have units of meters.

Curvilinear coordinates
~~~~~~~~~~~~~~~~~~~~~~~

Curvilinear coordinates are selected when the logical flag
**usingCurvilinearGrid** in namelist *PARM04* is set to true. The grid
spacing can not be set via the namelist. Instead, the grid descriptors
are read from data files, one for each descriptor. As for other grids,
the horizontal grid descriptors have units of meters.

Vertical grid
-------------

+----+----+
+----+----+

As for the horizontal grid, we use the suffixes “c” and “f” to indicates
faces and centers. Fig. [fig:vgrid]a shows the default vertical grid
used by the model. :math:`\Delta r_f` is the difference in :math:`r`
(vertical coordinate) between the faces (i.e. :math:`\Delta r_f \equiv -
\delta_k r` where the minus sign appears due to the convention that the
surface layer has index :math:`k=1`.).

The vertical grid is calculated in subroutine *INI\_VERTICAL\_GRID* and
specified via the vector **DELR** in namelist *PARM04*. The units of “r”
are either meters or Pascals depending on the isomorphism being used
which in turn is dependent only on the choice of equation of state.

There are alternative namelist vectors **DELZ** and **DELP** which
dictate whether z- or p- coordinates are to be used but we intend to
phase this out since they are redundant.

The reciprocals :math:`\Delta r_f^{-1}` and :math:`\Delta r_c^{-1}` are
pre-calculated (also in subroutine *INI\_VERTICAL\_GRID*). All vertical
grid descriptors are stored in common blocks in *GRID.h*.

The above grid (Fig. [fig:vgrid]a) is known as the cell centered
approach because the tracer points are at cell centers; the cell centers
are mid-way between the cell interfaces. This discretization is selected
when the thickness of the levels are provided (**delR**, parameter file
*data*, namelist *PARM04*) An alternative, the vertex or interface
centered approach, is shown in Fig. [fig:vgrid]b. Here, the interior
interfaces are positioned mid-way between the tracer nodes (no longer
cell centers). This approach is formally more accurate for evaluation of
hydrostatic pressure and vertical advection but historically the cell
centered approach has been used. An alternative form of subroutine
*INI\_VERTICAL\_GRID* is used to select the interface centered approach
This form requires to specify :math:`Nr+1` vertical distances **delRc**
(parameter file *data*, namelist *PARM04*, e.g.
*verification/ideal\_2D\_oce/input/data*) corresponding to surface to
center, :math:`Nr-1` center to center, and center to bottom distances.

Topography: partially filled cells
----------------------------------

<!– CMIREDIR:topo\_partial\_cells: –>

:raw-latex:`\cite{adcroft:97}` presented two alternatives to the
step-wise finite difference representation of topography. The method is
known to the engineering community as *intersecting boundary method*. It
involves allowing the boundary to intersect a grid of cells thereby
modifying the shape of those cells intersected. We suggested allowing
the topography to take on a piece-wise linear representation (shaved
cells) or a simpler piecewise constant representation (partial step).
Both show dramatic improvements in solution compared to the traditional
full step representation, the piece-wise linear being the best. However,
the storage requirements are excessive so the simpler piece-wise
constant or partial-step method is all that is currently supported.

Fig. [fig:hfacs] shows a schematic of the x-r plane indicating how the
thickness of a level is determined at tracer and u points. The physical
thickness of a tracer cell is given by :math:`h_c(i,j,k) \Delta
r_f(k)` and the physical thickness of the open side is given by
:math:`h_w(i,j,k) \Delta r_f(k)`. Three 3-D descriptors :math:`h_c`,
:math:`h_w` and :math:`h_s` are used to describe the geometry:
**hFacC**, **hFacW** and **hFacS** respectively. These are calculated in
subroutine *INI\_MASKS\_ETC* along with there reciprocals
**RECIP\_hFacC**, **RECIP\_hFacW** and **RECIP\_hFacS**.

The non-dimensional fractions (or h-facs as we call them) are calculated
from the model depth array and then processed to avoid tiny volumes. The
rule is that if a fraction is less than **hFacMin** then it is rounded
to the nearer of :math:`0` or **hFacMin** or if the physical thickness
is less than **hFacMinDr** then it is similarly rounded. The larger of
the two methods is used when there is a conflict. By setting
**hFacMinDr** equal to or larger than the thinnest nominal layers,
:math:`\min{(\Delta z_f)}`, but setting **hFacMin** to some small
fraction then the model will only lop thick layers but retain stability
based on the thinnest unlopped thickness;
:math:`\min{(\Delta z_f,\mbox{\bf hFacMinDr})}`.

Continuity and horizontal pressure gradient terms
=================================================

<!– CMIREDIR:continuity\_and\_horizontal\_pressure: –>

The core algorithm is based on the “C grid” discretization of the
continuity equation which can be summarized as:

.. math::

   \begin{aligned}
   \partial_t u + \frac{1}{\Delta x_c} \delta_i \left. \frac{ \partial \Phi}{\partial r}\right|_{s} \eta + \frac{\epsilon_{nh}}{\Delta x_c} \delta_i \Phi_{nh}' & = & G_u - \frac{1}{\Delta x_c} \delta_i \Phi_h' \label{eq:discrete-momu} \\
   \partial_t v + \frac{1}{\Delta y_c} \delta_j \left. \frac{ \partial \Phi}{\partial r}\right|_{s} \eta + \frac{\epsilon_{nh}}{\Delta y_c} \delta_j \Phi_{nh}' & = & G_v - \frac{1}{\Delta y_c} \delta_j \Phi_h' \label{eq:discrete-momv} \\
   \epsilon_{nh} \left( \partial_t w + \frac{1}{\Delta r_c} \delta_k \Phi_{nh}' \right) & = & \epsilon_{nh} G_w + \overline{b}^k - \frac{1}{\Delta r_c} \delta_k \Phi_{h}' \label{eq:discrete-momw} \\
   \delta_i \Delta y_g \Delta r_f h_w u +
   \delta_j \Delta x_g \Delta r_f h_s v +
   \delta_k {\cal A}_c w & = & {\cal A}_c \delta_k (P-E)_{r=0}
   \label{eq:discrete-continuity}\end{aligned}

 where the continuity equation has been most naturally discretized by
staggering the three components of velocity as shown in
Fig. [fig:cgrid3d]. The grid lengths :math:`\Delta x_c` and
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
   A}_c(P-E)_{r=0} \label{eq:discrete-freesurface}

 The source term :math:`P-E` on the rhs of continuity accounts for the
local addition of volume due to excess precipitation and run-off over
evaporation and only enters the top-level of the *ocean* model.

Hydrostatic balance
===================

<!– CMIREDIR:hydrostatic\_balance: –>

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
   \label{eq:discrete_hydro_ocean}

In the atmosphere, using p-coordinates, hydrostatic balance is
discretized:

.. math::

   \overline{\theta'}^k + \frac{1}{\Delta \Pi} \delta_k \Phi_h' = 0
   \label{eq:discrete_hydro_atmos}

 where :math:`\Delta \Pi` is the difference in Exner function between
the pressure points. The non-hydrostatic equations are not available in
the atmosphere.

The difference in approach between ocean and atmosphere occurs because
of the direct use of the ideal gas equation in forming the potential
energy conversion term :math:`\alpha \omega`. The form of these
conversion terms is discussed at length in
:raw-latex:`\cite{adcroft:02}`.

Because of the different representation of hydrostatic balance between
ocean and atmosphere there is no elegant way to represent both systems
using an arbitrary coordinate.

The integration for hydrostatic pressure is made in the positive
:math:`r` direction (increasing k-index). For the ocean, this is from
the free-surface down and for the atmosphere this is from the ground up.

The calculations are made in the subroutine *CALC\_PHI\_HYD*. Inside
this routine, one of other of the atmospheric/oceanic form is selected
based on the string variable **buoyancyRelation**.
