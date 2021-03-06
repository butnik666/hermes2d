#ifndef __HERMES2D_SEA_BREEZE_PARAMS_H
#define __HERMES2D_SEA_BREEZE_PARAMS_H

/* This file is generated by equations.py

   Just modify the python script if you want to change the quantities. Note
   that all the quantities in this file are related, so you *cannot* just
   change R or g (say) without also recalculating p, T and rho.
*/

// Everything is in SI units
// initial presure p in Pascals
#define p_z(z) (100000 - 11.3653418526412*(z))
// initial temperature T in Kelvin
#define T_z(z) (300.5 - 0.0341528522671867*(z))
// initial gas density
#define rho_z(z) (1.15894233531748)

// other physical constants
#define R 287.14           // Gas constant [J/(kg*K)]
#define g 0           // gravitational acceleration [m/s^2]
#define c_v 717.5       // specific heat capacity [J/(kg*K)]

// main characteristic constants
#define l_r 4000.0 // units: m
#define u_r 300.0 // units: m/s
#define rho_r 1.1 // units: kg/m^3
// other characteristic constants
#define t_r (l_r/u_r)  // time
#define p_r (rho_r*u_r*u_r) // presure
#define E_r (rho_r*u_r*u_r) // energy
#define g_r (l_r/(t_r*t_r)) // gravitational constant

// other constants
#define kappa (1 + R/c_v)

#endif
