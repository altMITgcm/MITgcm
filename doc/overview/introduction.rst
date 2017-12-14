Introduction
============

MITgcm has a number of novel aspects:

 - it can be used to study both atmospheric and oceanic phenomena; one hydrodynamical kernel is used to drive forward both atmospheric and oceanic models - see :numref:`onemodel`

  .. figure:: figs/onemodel.*
    :width: 80%
    :align: center
    :alt: One model for atmospheric and oceanic simulations
    :name: onemodel

    MITgcm has a single dynamical kernel that can drive forward either oceanic or atmospheric simulations.


 - it has a non-hydrostatic capability and so can be used to study both small-scale and large scale processes - see :numref:`all-scales`

  .. figure:: figs/scales.*
    :width: 90%
    :align: center
    :alt: MITgcm can simulate a wide range of scales
    :name: all-scales

    MITgcm has non-hydrostatic capabilities, allowing the model to address a wide range of phenomenon - from convection on the left, all the way through to global circulation patterns on the right.

 - finite volume techniques are employed yielding an intuitive discretization and support for the treatment of irregular geometries using orthogonal curvilinear grids and shaved cells - see :numref:`fvol`

  .. figure:: figs/fvol.*
    :width: 80%
    :align: center
    :alt: Finit volume techniques
    :name: fvol

    Finite volume techniques (bottom panel) are used, permitting a treatment of topography that rivals :math:`\sigma` (terrain following) coordinates.

 - tangent linear and adjoint counterparts are automatically maintained along with the forward model, permitting sensitivity and optimization studies.

 - the model is developed to perform efficiently on a wide variety of computational platforms.


Key publications reporting on and charting the development of the model are Hill and Marshall (1995), Marshall et al. (1997a), 
Marshall et al. (1997b), Adcroft and Marshall (1997), Marshall et al. (1998), Adcroft and Marshall (1999), Hill et al. (1999),
Marotzke et al. (1999), Adcroft and Campin (2004), Adcroft et al. (2004b), Marshall et al. (2004) (an overview on the model formulation can also be found in Adcroft et al. (2004c)):

Hill, C. and J. Marshall, (1995)
Application of a Parallel Navier-Stokes Model to Ocean Circulation in 
Parallel Computational Fluid Dynamics,
In Proceedings of Parallel Computational Fluid Dynamics: Implementations 
and Results Using Parallel Computers, 545-552.
Elsevier Science B.V.: New York :cite:`hill:95`

Marshall, J., C. Hill, L. Perelman, and A. Adcroft, (1997a)
Hydrostatic, quasi-hydrostatic, and nonhydrostatic ocean modeling,
J. Geophysical Res., **102(C3)**, 5733-5752 :cite:`marshall:97a`

Marshall, J., A. Adcroft, C. Hill, L. Perelman, and C. Heisey, (1997b)
A finite-volume, incompressible Navier Stokes model for studies of the ocean
on parallel computers, J. Geophysical Res., **102(C3)**, 5753-5766 :cite:`marshall:97b`

Adcroft, A.J., Hill, C.N. and J. Marshall, (1997)
Representation of topography by shaved cells in a height coordinate ocean
model, Mon Wea Rev, **125**, 2293-2315 :cite:`adcroft:97`

Marshall, J., Jones, H. and C. Hill, (1998)
Efficient ocean modeling using non-hydrostatic algorithms,
Journal of Marine Systems, **18**, 115-134 :cite:`mars-eta:98`

Adcroft, A., Hill C. and J. Marshall: (1999)
A new treatment of the Coriolis terms in C-grid models at both high and low
resolutions,
Mon. Wea. Rev., **127**, 1928-1936 :cite:`adcroft:99`

Hill, C, Adcroft,A., Jamous,D., and J. Marshall, (1999)
A Strategy for Terascale Climate Modeling,
In Proceedings of the Eighth ECMWF Workshop on the Use of Parallel Processors
in Meteorology, 406-425
World Scientific Publishing Co: UK :cite:`hill:99`

Marotzke, J, Giering,R., Zhang, K.Q., Stammer,D., Hill,C., and T.Lee, (1999)
Construction of the adjoint MIT ocean general circulation model and 
application to Atlantic heat transport variability,
J. Geophysical Res., **104(C12)**, 29,529-29,547 :cite:`maro-eta:99`

A. Adcroft and J.-M. Campin, (2004a)
Re-scaled height coordinates for accurate representation of free-surface flows in ocean circulation models, 
Ocean Modelling, **7**, 269–284 :cite:`adcroft:04a`

A. Adcroft, J.-M. Campin, C. Hill, and J. Marshall, (2004b)
Implementation of an atmosphere-ocean general circulation model on the expanded 
spherical cube, 
Mon Wea Rev , **132**, 2845–2863 :cite:`adcroft:04b`

J. Marshall, A. Adcroft, J.-M. Campin, C. Hill, and A. White, (2004)
Atmosphere-ocean modeling exploiting fluid isomorphisms, Mon. Wea. Rev., **132**, 2882–2894 :cite:`marshall:04`

A. Adcroft, C. Hill, J.-M. Campin, J. Marshall, and P. Heimbach, (2004c)
Overview of the formulation and numerics of the MITgcm, In Proceedings of the ECMWF seminar series on Numerical Methods, Recent developments in numerical methods for atmosphere and ocean modelling, 139–149. URL: http://mitgcm.org/pdfs/ECMWF2004-Adcroft.pdf :cite:`adcroft:04c`

We begin by briefly showing some of the results of the model in action to
give a feel for the wide range of problems that can be addressed using it.
