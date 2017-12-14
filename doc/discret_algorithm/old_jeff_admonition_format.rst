.. admonition:: Pressure method calling tree
  :class: note

  - :filelink:`FORWARD\_STEP <model/src/forward_step.F>`
    
    - :filelink:`DYNAMICS <model/src/dynamics.F>`
      
      - :filelink:`TIMESTEP <model/src/timestep.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxxx}` :math:`u^*,v^*` :eq:`ustar-rigid-lid` , :eq:`vstar-rigid-lid`
        
    - :filelink:`SOLVE\_FOR\_PRESSURE <model/src/solve_for_pressure.F>`
         
      - :filelink:`CALC\_DIV\_GHAT <model/src/calc_div_ghat.F>` :math:`\phantom{xxxxxxxxxxxxxxxx}` :math:`H\widehat{u^*},H\widehat{v^*}` :eq:`elliptic`
           
      - :filelink:`CG2D <model/src/cg2d.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxxxxxxxxx}` :math:`\eta^{n+1}` :eq:`elliptic`
       
    - :filelink:`MOMENTUM\_CORRECTION\_STEP <model/src/momentum_correction_step.F>`
        
      - :filelink:`CALC\_GRAD\_PHI\_SURF <model/src/calc_grad_phi_surf.F>` :math:`\phantom{xxxxxxxxxx}` :math:`\nabla \eta^{n+1}`
        
      - :filelink:`CORRECTION\_STEP  <model/src/correction_step.F>` :math:`\phantom{xxxxxxxxxxxxw}` :math:`u^{n+1},v^{n+1}` :eq:`un+1-rigid-lid` , :eq:`vn+1-rigid-lid`

.. admonition:: Adams-Bashforth calling tree
  :class: note

  - :filelink:`FORWARD\_STEP <model/src/forward_step.F>`
    
    - :filelink:`THERMODYNAMICS <model/src/thermodynamics.F>`
      
      - :filelink:`TEMP\_INTEGRATE <model/src/temp_integrate.F>`
        
        - :filelink:`GAD\_CALC\_RHS <pkg/generic_advdiff/gad_calc_rhs.F>` :math:`\phantom{xxxxxxxxxw}` :math:`G_\theta^n = G_\theta( u, \theta^n)`
          
        - either 

          - :filelink:`EXTERNAL\_FORCING <model/src/external_forcing.F>` :math:`\phantom{xxx}` :math:`G_\theta^n = G_\theta^n + {\cal Q}`
            
          - :filelink:`ADAMS\_BASHFORTH2 <model/src/adams_bashforth2.F>` :math:`\phantom{xi}` :math:`G_\theta^{(n+1/2)}` :eq:`adams-bashforth2`
        
        - or

          - :filelink:`EXTERNAL\_FORCING <model/src/external_forcing.F>` :math:`\phantom{xxx}` :math:`G_\theta^{(n+1/2)} = G_\theta^{(n+1/2)} + {\cal Q}`
        
      - :filelink:`TIMESTEP\_TRACER <model/src/timestep_tracer.F>` :math:`\phantom{xxxxxxxxxxx}` :math:`\tau^*` :eq:`taustar`
        
      - :filelink:`IMPLDIFF  <model/src/impldiff.F>` :math:`\phantom{xxxxxxxxxxxxxxxxxxw}` :math:`\tau^{(n+1)}` :eq:`tau-n+1-implicit`


