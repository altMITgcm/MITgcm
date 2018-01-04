
Horizontal grid
---------------

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

  .. figure:: figs/hgrid-abcd.*
    :width: 90%
    :align: center
    :alt: hgrid-abcd
    :name: hgrid-abcd
       Staggering of horizontal grid descriptors (lengths and areas). The grid lines indicate the tracer cell boundaries and are the reference grid for all panels. a) The area of a tracer cell, :math:`A_c`, is bordered by the lengths :math:`\Delta x_g` and :math:`\Delta y_g`. b) The area of a vorticity cell, :math:`A_\zeta`, is bordered by the lengths :math:`\Delta x_c` and :math:`\Delta y_c`. c) The area of a u cell, :math:`A_w`, is bordered by the lengths :math:`\Delta x_v` and :math:`\Delta y_f`. d) The area of a v cell, :math:`A_s`, is bordered by the lengths :math:`\Delta x_f` and :math:`\Delta y_u`.

:numref:`hgrid-abcd` (a) shows the tracer cell (synonymous with the continuity
cell). The length of the southern edge, :math:`\Delta x_g`, western
edge, :math:`\Delta y_g` and surface area, :math:`A_c`, presented in the
vertical are stored in arrays **DXg**, **DYg** and **rAc**. The “g”
suffix indicates that the lengths are along the defining grid
boundaries. The “c” suffix associates the quantity with the cell
centers. The quantities are staggered in space and the indexing is such
that **DXg(i,j)** is positioned to the south of **rAc(i,j)** and
**DYg(i,j)** positioned to the west.

:numref:`hgrid-abcd` (b) shows the vorticity cell. The length of the southern
edge, :math:`\Delta x_c`, western edge, :math:`\Delta y_c` and surface
area, :math:`A_\zeta`, presented in the vertical are stored in arrays
**DXc**, **DYc** and **rAz**. The “z” suffix indicates that the lengths
are measured between the cell centers and the “:math:`\zeta`” suffix
associates points with the vorticity points. The quantities are
staggered in space and the indexing is such that **DXc(i,j)** is
positioned to the north of **rAz(i,j)** and **DYc(i,j)** positioned to
the east.

:numref:`hgrid-abcd` (c) shows the “u” or western (w) cell. The length of the
southern edge, :math:`\Delta x_v`, eastern edge, :math:`\Delta y_f` and
surface area, :math:`A_w`, presented in the vertical are stored in
arrays **DXv**, **DYf** and **rAw**. The “v” suffix indicates that the
length is measured between the v-points, the “f” suffix indicates that
the length is measured between the (tracer) cell faces and the “w”
suffix associates points with the u-points (w stands for west). The
quantities are staggered in space and the indexing is such that
**DXv(i,j)** is positioned to the south of **rAw(i,j)** and **DYf(i,j)**
positioned to the east.

:numref:`hgrid-abcd` (d) shows the “v” or southern (s) cell. The length of the
northern edge, :math:`\Delta x_f`, western edge, :math:`\Delta y_u` and
surface area, :math:`A_s`, presented in the vertical are stored in
arrays **DXf**, **DYu** and **rAs**. The “u” suffix indicates that the
length is measured between the u-points, the “f” suffix indicates that
the length is measured between the (tracer) cell faces and the “s”
suffix associates points with the v-points (s stands for south). The
quantities are staggered in space and the indexing is such that
**DXf(i,j)** is positioned to the north of **rAs(i,j)** and **DYu(i,j)**
positioned to the west.

.. admonition:: S/R :filelink:`INI_CARTESIAN_GRID <model/src/ini_cartesian_grid.F>` , :filelink:`INI_SPHERICAL_POLAR_GRID <model/src/ini_spherical_polar_grid.F>` , :filelink:`INI_CURVILINEAR_GRID <model/src/ini_curvilinear_grid.F>`
  :class: note

    | :math:`A_c , A_\zeta , A_w , A_s :` : **rAc, rAz, rAw, rAs** ( :filelink:`GRID.h <model/inc/GRID.h>` )
    | :math:`\Delta x_g , \Delta y_g` : **DXg, DYg** ( :filelink:`GRID.h <model/inc/GRID.h>` )
    | :math:`\Delta x_c , \Delta y_c` : **DXc, DYc** ( :filelink:`GRID.h <model/inc/GRID.h>` )
    | :math:`\Delta x_f , \Delta y_f` : **DXf, DYf** ( :filelink:`GRID.h <model/inc/GRID.h>` )
    | :math:`\Delta x_v , \Delta y_u` : **DXv, DYu** ( :filelink:`GRID.h <model/inc/GRID.h>` )


Reciprocals of horizontal grid descriptors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lengths and areas appear in the denominator of expressions as much as in
the numerator. For efficiency and portability, we pre-calculate the
reciprocal of the horizontal grid quantities so that in-line divisions
can be avoided.

For each grid descriptor (array) there is a reciprocal named using the
prefix ``RECIP_``. This doubles the amount of storage in :filelink:`GRID.h <model/inc/GRID.h>` but
they are all only 2-D descriptors.

.. admonition:: S/R :filelink:`INI_MASKS_ETC <model/src/ini_masks_etc.F>`
  :class: note

    | :math:`A_c^{-1} , A_\zeta^{-1} , A_w^{-1} , A_s^{-1} :` : **RECIP_Ac, RECIP_Az, RECIP_Aw, RECIP_As** ( :filelink:`GRID.h <model/inc/GRID.h>` )
    | :math:`\Delta x_g^{-1} , \Delta y_g^{-1}` : **RECIP_DXg, RECIP_DYg** ( :filelink:`GRID.h <model/inc/GRID.h>` )
    | :math:`\Delta x_c^{-1} , \Delta y_c^{-1}` : **RECIP_DXc, RECIP_DYc** ( :filelink:`GRID.h <model/inc/GRID.h>` )
    | :math:`\Delta x_f^{-1} , \Delta y_f^{-1}` : **RECIP_DXf, RECIP_DYf** ( :filelink:`GRID.h <model/inc/GRID.h>` )
    | :math:`\Delta x_v^{-1} , \Delta y_u^{-1}` : **RECIP_DXv, RECIP_DYu** ( :filelink:`GRID.h <model/inc/GRID.h>` )

Cartesian coordinates
~~~~~~~~~~~~~~~~~~~~~

Cartesian coordinates are selected when the logical flag
:varlink:`usingCartesianGrid` in namelist ``PARM04`` is set to true. The grid
spacing can be set to uniform via scalars :varlink:`dXspacing` and
:varlink:`dYspacing` in namelist ``PARM04`` or to variable resolution by the
vectors :varlink:`DELX` and :varlink:`DELY`. Units are normally meters.
Non-dimensional coordinates can be used by interpreting the
gravitational constant as the Rayleigh number.

Spherical-polar coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spherical coordinates are selected when the logical flag
:varlink:`usingSphericalPolarGrid` in namelist ``PARM04`` is set to true. The
grid spacing can be set to uniform via scalars :varlink:`dXspacing` and
:varlink:`dYspacing` in namelist ``PARM04`` or to variable resolution by the
vectors :varlink:`DELX` and :varlink:`DELY`. Units of these namelist variables are
alway degrees. The horizontal grid descriptors are calculated from these
namelist variables have units of meters.

Curvilinear coordinates
~~~~~~~~~~~~~~~~~~~~~~~

Curvilinear coordinates are selected when the logical flag
:varlink:`usingCurvilinearGrid` in namelist ``PARM04`` is set to true. The grid
spacing can not be set via the namelist. Instead, the grid descriptors
are read from data files, one for each descriptor. As for other grids,
the horizontal grid descriptors have units of meters.

