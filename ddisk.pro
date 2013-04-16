;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ddisk.pro
;
;Version 4.0, 5-March-2010.
;
;A circular planetesimal rings generate dust along discreet streamlines distributed uniformly
;about the ring. The resulting dusty orbits are eccentric due to radiation pressure, and are
;determined by their beta parameter = radiation pressure/stellar gravity. The
;abundance of dust n(t) in each orbit evolves over time due to dust production and collisions
;with other dust grains. This code solve for the number of grains in each in orbits and over
;time, assuming a power-law for the ring's dust production rate versus grain size.

;By Joseph M. Hahn, Space Science Institute, jhahn@spacescience.org.
;Feel free to email me any suggestions or bug reports. See 'Diagnosing Circumstellar Debris
;Disks' (Hahn, 2010, ApJ) for more information.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Clock time (in seconds) when execution starts
t_start = systime( /seconds)

;Compile the subroutines that are called below
.run odeint, mpfitfun, mpfit, subroutines

;Read input parameters.
@inputs

;Define various constants, and store them in structure constants. Values from
;Astrophysical Data (Lang, 1992).
L_sun_cgs     = 3.8268e33                      ;solar luminosity in ergs/sec
G_cgs         = 6.6726e-8                      ;gravitational constant in cm^3/gm/sec^2
M_sun_gm      = 1.9891e33                      ;solar mass in gm
M_earth_gm    = 5.9742e27                      ;Earth mass in gm
M_moon_gm     = 7.3483e25                      ;Moon mass in gm
c_cgs         = 2.9979e10                      ;speed of light in cm/sec
h_cgs         = 6.6261d-27                     ;Planck's constant in ergs*sec
k_cgs         = 1.3807d-16                     ;Boltzmann constant in erg/K
sbc_cgs       = 5.6705d-5                      ;Stefan-Boltzmann cnst, erg/cm^2/K^4/sec
Jansky        = 1.0000e-23                     ;one Jansky in units of ergs/cm^2/sec/Hz
cm_per_au     = 1.4960e13                      ;number of cm in 1 AU
km_per_au     = 1.4960e8                       ;number of km in 1 AU
nm_per_cm     = 1.0000d7                       ;number of nanometers in one cm 
micron_per_cm = 1.0000e4                       ;number of microns in a cm
sec_per_yr    = 365.25*24.0*60.0*60.0          ;seconds in one yr
twopi         = 2d*!dpi                        ;2*Pi
deg2_per_ster = (360.0/twopi)^2                ;number of degrees in one steradian
arc2_per_ster = (3600d^2)*deg2_per_ster        ;number of arcsec^2 in one steradian
constants = { L_sun_cgs:L_sun_cgs, G_cgs:G_cgs, M_sun_gm:M_sun_gm, M_earth_gm:M_earth_gm, $
    M_moon_gm:M_moon_gm, c_cgs:c_cgs, h_cgs:h_cgs, k_cgs:k_cgs, sbc_cgs:sbc_cgs, $
    Jansky:Jansky, cm_per_au:cm_per_au, km_per_au:km_per_au, nm_per_cm:nm_per_cm, $
    micron_per_cm:micron_per_cm, sec_per_yr:sec_per_yr, twopi:twopi, $
    deg2_per_ster:deg2_per_ster, arc2_per_ster:arc2_per_ster }

;Define psym=8 as filled dots
Nd = 41
phi_d = twopi*findgen(Nd)/(Nd - 1)
usersym, cos(phi_d), sin(phi_d), fill=1

;Generate beta, orbit elements, and dust sizes for all streamlines
make_orbits, inputs, ar, a_ring, ring_idx, beta, beta_idx, a, e, Q_peri, Q_apo

;Calculate alpha_star=scaled collisional probability density
calc_alpha, inputs, constants, a, e, beta, ar, a_ring, ring_idx, alpha_star, $
    r_crossing, N_crossing, N_collisions, N_disruptions, grazing_orbits=0
t_stop = systime( /seconds)
execution_time_min = (t_stop - t_start)/60.0
mx = max(alpha_star, idx)
idx = array_indices(alpha_star, idx)
print, 'max alpha, beta_target, beta_impactor          = ', mx, beta(idx) 
print, 'execution time (minutes)                       = ', execution_time_min

;Solve for the grains' equilibrium populations n_star_eq
equilibrium, alpha_star, n_star_eq 

;Evolve collision evolution equations over dimensionless time tau (which is called t^star
;in the paper), if desired
n_star = 0
tau = 0.0
if (evolve_system eq 'YES') then $
    evolve_system, alpha_star, N_times, N_longitudes, tau1, tau2, tau_stop, n_star, tau

;Save results
t_stop = systime( /seconds)
execution_time_min = (t_stop - t_start)/60.0
print, 'execution time (minutes)                       = ', execution_time_min
save, filename='ddisk.dat', alpha_star, ar, n_star, n_star_eq, tau, $
    a, e, beta, beta_idx, a_ring, ring_idx, Q_apo, Q_peri, execution_time_min, $
    r_crossing, N_crossing, constants, inputs, N_collisions, N_disruptions
