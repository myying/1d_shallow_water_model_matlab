1D Shallow Water Equation  - Yue (Michael) Ying 2015
For the purpose of testing basic numerical schemes for advection-diffusion problems.

Recoded to MATLAB from a Fortran version provided by Meteo 526 NWP class.
Setting the following parameters as in "param.ctr" from the Fortran code:
  C (same as "vel")
  hbar
  amp
  dx
  dt
  g
  f
  kdif (same as "K"), 
  hc (same as "(jhin+jhend)/2")
  hw (same as "jhend-jhin")
  nifcor
  nifwind
  nifdif
  nifad
  
Using "advection_scheme=1" will do the same thing as the Fortran code, other options 
are experimental.
