Shapiro Filter
==============

The Shapiro filter :raw-latex:`\cite{Shapiro_70}` is a high order
horizontal filter that efficiently remove small scale grid noise without
affecting the physical structures of a field. It is applied at the end
of the time step on both velocity and tracer fields.

Three different space operators are considered here (S1,S2 and S4). They
differs essentially by the sequence of derivative in both X and Y
directions. Consequently they show different damping response function
specially in the diagonal directions X+Y and X-Y.

Space derivatives can be computed in the real space, taken into account
the grid spacing. Alternatively, a pure computational filter can be
defined, using pure numerical differences and ignoring grid spacing.
This later form is stable whatever the grid is, and therefore specially
useful for highly anisotropic grid such as spherical coordinate grid. A
damping time-scale parameter :math:`\tau_{shap}` defines the strength of
the filter damping.

The 3 computational filter operators are :

.. math::

   \mathrm{S1c:}\hspace{2cm}
   [1 - 1/2 \frac{\Delta t}{\tau_{shap}}
      \{ (\frac{1}{4}\delta_{ii})^n 
       + (\frac{1}{4}\delta_{jj})^n \} ]

.. math::

   \mathrm{S2c:}\hspace{2cm}
   [1 - \frac{\Delta t}{\tau_{shap}} 
   \{ \frac{1}{8} (\delta_{ii} + \delta_{jj}) \}^n]

.. math::

   \mathrm{S4c:}\hspace{2cm}
   [1 - \frac{\Delta t}{\tau_{shap}} (\frac{1}{4}\delta_{ii})^n]
   [1 - \frac{\Delta t}{\tau_{shap}} (\frac{1}{4}\delta_{jj})^n]

In addition, the S2 operator can easily be extended to a physical space
filter:

.. math::

   \mathrm{S2g:}\hspace{2cm}
   [1 - \frac{\Delta t}{\tau_{shap}} 
   \{ \frac{L_{shap}^2}{8} \overline{\nabla}^2 \}^n]

with the Laplacian operator :math:`\overline{\nabla}^2 ` and a length
scale parameter :math:`L_{shap}`. The stability of this S2g filter
requires :math:`L_{shap} < \mathrm{Min}^{(Global)}(\Delta x,\Delta y)`.

SHAP Diagnostics
----------------

::


    ------------------------------------------------------------------------
    <-Name->|Levs|<-parsing code->|<--  Units   -->|<- Tile (max=80c) 
    ------------------------------------------------------------------------
    SHAP_dT |  5 |SM      MR      |K/s             |Temperature Tendency due to Shapiro Filter
    SHAP_dS |  5 |SM      MR      |g/kg/s          |Specific Humidity Tendency due to Shapiro Filter
    SHAP_dU |  5 |UU   148MR      |m/s^2           |Zonal Wind Tendency due to Shapiro Filter
    SHAP_dV |  5 |VV   147MR      |m/s^2           |Meridional Wind Tendency due to Shapiro Filter
