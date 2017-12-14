Notation
========

Because of the particularity of the vertical direction in stratified
fluid context, in this chapter, the vector notations are mostly used for
the horizontal component: the horizontal part of a vector is simply
written :math:`\vec{\bf v}` (instead of :math:`{\bf v_h}` or
:math:`\vec{\mathbf{v}}_{h}` in chaper 1) and a 3.D vector is simply
written :math:`\vec{v}` (instead of :math:`\vec{\mathbf{v}}` in chapter
1).

| The notations we use to describe the discrete formulation of the model
  are summarized hereafter:
| general notation:
| :math:`\Delta x, \Delta y, \Delta r` grid spacing in X,Y,R directions.
| :math:`A_c,A_w,A_s,A_{\zeta}` : horizontal area of a grid cell
  surrounding :math:`\theta,u,v,\zeta` point.
| :math:`{\cal V}_u , {\cal V}_v , {\cal V}_w , {\cal V}_\theta` :
  Volume of the grid box surrounding :math:`u,v,w,\theta` point;
| :math:`i,j,k` : current index relative to X,Y,R directions;
| basic operator:
| :math:`\delta_i ` :
  :math:`\delta_i \Phi = \Phi_{i+1/2} - \Phi_{i-1/2} ` [eq:delta\_i]
| :math:`~^{-i}` :
  :math:`\overline{\Phi}^i = ( \Phi_{i+1/2} + \Phi_{i-1/2} ) / 2 `
  [eq:bar\_i]
| :math:`\delta_x ` :
  :math:`\delta_x \Phi = \frac{1}{\Delta x} \delta_i \Phi `
  [eq:delta\_x]
| :math:`\overline{\nabla}` = horizontal gradient operator :
  :math:`\overline{\nabla} \Phi = \{ \delta_x \Phi , \delta_y \Phi \}`
  [eq:d\_grad]
| :math:`\overline{\nabla} \cdot` = horizontal divergence operator :
  :math:`\overline{\nabla}\cdot \vec{\mathrm{f}}  = 
  \frac{1}{\cal A} \{ \delta_i \Delta y \, \mathrm{f}_x 
                    + \delta_j \Delta x \, \mathrm{f}_y \} ` [eq:d\_div]
| :math:`\overline{\nabla}^2 ` = horizontal Laplacian operator :
  :math:` \overline{\nabla}^2 \Phi = 
     \overline{\nabla}\cdot \overline{\nabla}\Phi ` [eq:d\_lap]
