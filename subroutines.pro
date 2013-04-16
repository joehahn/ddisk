;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;subroutines.pro
;
;Procedures called by v4.0 of ddisk.pro,
;by Joseph M. Hahn, Space Science Institute, jhahn@spacescience.org, 5 March 2010.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Execute setplot to direct plots to a postscript file. Then make your figures,
;and then run output_plot.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setplot, large=large, lr=lr, narrow=narrow, short=short, two_by_two=two_by_two

set_plot,'ps'
device,xsize=16,ysize=15.2,xoffset=3,yoffset=8,/portrait,bits_per_pixel=8,/color
if (keyword_set(large) eq 1) then $
  device,xsize=17.5,ysize=24,xoffset=2,yoffset=2,/color,/portrait,bits_per_pixel=8
if (keyword_set(lr) eq 1) then $
  device,xsize=24,ysize=10.8,xoffset=2,yoffset=6,/color,/portrait,bits_per_pixel=8
if (keyword_set(narrow) eq 1) then $
  device,xsize=16,ysize=8,xoffset=3,yoffset=8,/portrait,bits_per_pixel=8,/color
if (keyword_set(short) eq 1) then $
  device,xsize=20.0,ysize=10.0,xoffset=2,yoffset=2,/color,/portrait,bits_per_pixel=8
if (keyword_set(two_by_two) eq 1) then $
  device,xsize=20.0,ysize=18.2,xoffset=2,yoffset=4,/color,/portrait,bits_per_pixel=8

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;This moves file 'idl.ps' to desired file 'newfile', and directs plots back to IDL's
;window on the monitor.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro output_plot, new_file

device,/close
file='idl.'+strlowcase(!d.name)
if n_elements(new_file) eq 0 then spawn,'lpr '+file
if n_elements(new_file) eq 0 then spawn,'rm idl.ps' $
else begin
    spawn,'mv '+file+' '+new_file
    print,'image sent to file ',new_file
endelse
set_plot,'x'
return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Calculate the polar coordinates (r_cross, theta_cross) where two coplanar ellipses
;cross. Two solutions are usually possible, so the outputs r_cross and theta_cross are
;2-element arrays. If orbits don't cross, then r_cross[k] = -1 is returned. Inputs
;are the two orbits semimajor axes a1,a2, eccentricities e1,e2, and longitudes of
;periapse w1,w2 (in radians).
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro orbit_crossing, a1, e1, w1, a2, e2, w2, r_c, theta_c

;Storage for crossing sites. If crossing site is not found, then r_c[j]=-1 
r_c = [-1.0, -1.0]
theta_c = [0.0, 0.0]

;Semilatus rectum
p1 = a1*(1 - e1^2)
p2 = a2*(1 - e2^2)
p = p1/p2

;relative longitude of periapse
w = w2 - w1

;various constants
a = e1 - p*e2*cos(w)
b = p*e2*sin(w)
c = 1 - p
a2_plus_b2  = a^2 + b^2
;Return if a=0=b, since these orbits don't cross
if ((a eq 0) and (b eq 0)) then return
BB = a*c/a2_plus_b2
CC = (c^2 - b^2)/a2_plus_b2
BB2_minus_CC = BB^2 - CC

;A small single-precision floating point number that is larger than roundoff error.
;All angles that differ by less than this amount (in radians) are considered equal.
small = 1.0e-5

;Special case for when b=0 or BB2_minus_CC=0 to within a small numerical error
if ( abs(b) le small ) then begin $
    cos_theta = -BB & $
    ;If |cos_theta| is a smidge above 1 due to roundoff error, set |cos_theta|=1-small/2
    diff = abs(cos_theta) - 1 & $
    if (abs(diff) le small) then begin $
        if (cos_theta gt 0) then cos_theta =  1.0 - 0.5*small & $
        if (cos_theta lt 0) then cos_theta = -1.0 + 0.5*small & $
    endif & $
    ;orbits cross if |cos(theta)| <= 1
    if (abs(cos_theta) le 1) then begin $
        theta = acos(cos_theta) & $
        theta_c = [ theta, -theta ] & $
        r_c = p1/(1 + e1*cos(theta_c)) & $
        theta_c = theta_c + w1 & $
    endif & $
endif else begin $

    ;Otherwise orbits cross if BB2_minus_CC>0, |sin(theta)|<=1, and |cos(theta)|<=1
    if ( BB2_minus_CC ge 0) then begin $

        ;Obtain first solution, and store if valid
        cos_theta = -BB + sqrt(BB2_minus_CC) & $
        sin_theta = (a*cos_theta + c)/b & $
        ;If |cos_theta| or |sin_theta| is a smidge above 1 due to roundoff,
        ;set to 1 - small/2
        diff = abs(cos_theta) - 1 & $
        if (abs(diff) le small) then begin $
            if (cos_theta gt 0) then cos_theta =  1.0 - 0.5*small & $
            if (cos_theta lt 0) then cos_theta = -1.0 + 0.5*small & $
        endif & $
        diff = abs(sin_theta) - 1 & $
        if (abs(diff) le small) then begin $
            if (sin_theta gt 0) then sin_theta =  1.0 - 0.5*small & $
            if (sin_theta lt 0) then sin_theta = -1.0 + 0.5*small & $
        endif & $
        if ( (abs(cos_theta) le 1) and (abs(sin_theta) le 1) ) then begin $
            theta = atan(sin_theta, cos_theta) & $
            r = p1/(1 + e1*cos(theta)) & $
            theta = theta + w1 & $
            r_c[0] = r & $
            theta_c[0] = theta & $
        endif & $

        ;Obtain second solution, and store if valid
        cos_theta = -BB - sqrt(BB2_minus_CC) & $
        sin_theta = (a*cos_theta + c)/b & $
        ;If |cos_theta| or |sin_theta| is a smidge above 1 due to roundoff,
        ;set to 1 - small/2
        diff = abs(cos_theta) - 1 & $
        if (abs(diff) le small) then begin $
            if (cos_theta gt 0) then cos_theta =  1.0 - 0.5*small & $
            if (cos_theta lt 0) then cos_theta = -1.0 + 0.5*small & $
        endif & $
        diff = abs(sin_theta) - 1 & $
        if (abs(diff) le small) then begin $
            if (sin_theta gt 0) then sin_theta =  1.0 - 0.5*small & $
            if (sin_theta lt 0) then sin_theta = -1.0 + 0.5*small & $
        endif & $
        if ( (abs(cos_theta) le 1) and (abs(sin_theta) le 1) ) then begin $
            theta = atan(sin_theta, cos_theta) & $
            r = p1/(1 + e1*cos(theta)) & $
            theta = theta + w1 & $
            r_c[1] = r & $
            theta_c[1] = theta & $
        endif & $

    endif & $
endelse

;If first solution = second solution to within numerical error, then the two
;orbits are grazing, so eliminate the redundant solution.
diff = sqrt( (r_c[0] - r_c[1])^2 + (theta_c[0] - theta_c[1])^2 )
if (diff le small) then begin $
    r_c[1] = -1.0 & $
    theta_c[1] = 0.0 & $
endif

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Calculate radial (v_r), tangential (v_theta), vertical (v_z), and total (v) velocity
;of grain having orbit elements a,e,w and size beta at polar coordinates (r, theta), in
;units of the outermost planetesimal ring velocity Vr. Note that ar=1=ring radius,
;and that only a mean vertical velocity v_z=I*v is calculated, 
;where I=dust inclination in radians.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro velocity, a, e, w, beta, r, theta, ar, I, v_r, v_theta, v_z, v

na = sqrt((ar/a)*(1.0 - beta))
factor = na/sqrt(1.0 - e^2)
anomaly = theta - w
v_r = factor*e*sin(anomaly)
v_theta = factor*(1.0 + e*cos(anomaly))
v = na*sqrt(2.0*a/r - 1.0)
v_z = I*v

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Make orbits. Generate orbit elements a,e and peri/apo distances Q_peri, Q_apo. All
;a's are in units of the outermost planetesimal ring's semimajor axis, which is alway ar=1.
;Also generate the planetesimal rings' semimajor axis a_ring. Use ring_idx to indicate which
;planetesimal ring manufactures which dust orbit such that ring_idx[5]=2 means that ring #2
;is responsible for making the dust in orbit #5.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro make_orbits, inputs, ar, a_ring, ring_idx, beta, beta_idx, a, e, Q_peri, Q_apo

;Extract inputs
Q_apo_min = inputs.Q_apo_min
Q_apo_max = inputs.Q_apo_max
N_sizes = inputs.N_sizes
N_rings = inputs.N_rings
ar_inner = inputs.ar_inner

;The planetesimal disk's outer radius ar is the system's scale length,
;and is always unity.
ar = 1.0

;Minimum/maximum beta for dust generated by the planetesimal rings.
beta_min = 0.5*(1 - ar/Q_apo_min)
beta_max = 0.5*(1 - ar/Q_apo_max)

;Dust radii in units of s_min
s_min = 1.0
s_max = beta_max/beta_min
s0 = (s_max - s_min)*findgen(N_sizes)/(N_sizes - 1) + s_min
s0 = rotate(s0, 2)

;Dimensionless size parameter beta.
beta0 = beta_max/s0

;Planetesimal ring semimajor axes
a_ring0 = 1.0
if (N_rings gt 1) then $
    a_ring0 = (ar - ar_inner)*findgen(N_rings)/(N_rings - 1) + ar_inner

;Integers are also used to identify grains by their size and their ring of origin
beta0_idx = indgen(N_sizes)
ring0_idx = indgen(N_rings)

;Create a 2D array for beta[size-index, ring-index]. Ditto for beta_idx and a_ring
beta = fltarr(N_sizes, N_rings)
beta_idx = indgen(N_sizes, N_rings)
a_ring = fltarr(N_sizes, N_rings)
ring_idx = indgen(N_sizes, N_rings)
for j = 0, N_rings - 1 do begin $
    beta[0, j] = beta0 & $
    beta_idx[0, j] = beta0_idx & $
    a_ring[*, j] = a_ring0[j] & $
    ring_idx[*, j] = ring0_idx[j] & $
endfor

;Convert to 1D arrays
beta = beta[*]
beta_idx = beta_idx[*]
a_ring = a_ring[*]
ring_idx = ring_idx[*]

;Orbit elements, periapse, and apoapse
a = (1 - beta)*a_ring/(1 - 2*beta)
e = beta/(1 - beta)
Q_peri = a*(1 - e)
Q_apo  = a*(1 + e)

;Plot streamlines periapse and apoapose distance versus beta
plot_io, beta, Q_apo, psym=8, yrange=[0.8*min(Q_peri), 1.2*max(Q_apo)], $
    xtitle='!7b!3', charsize=1.5, ytitle='q!lperi!n and Q!lapo!n    (a!lr!n)', $
    symsize=0.8, ystyle=1 & $
for j = 0, N_rings - 1 do begin $
    k = where(ring_idx eq j) & $
    oplot, beta[k], Q_peri[k], color=128 & $
    oplot, beta[k], Q_apo[k], color=128 & $
endfor

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Calculate scaled collisional probability density alpha_star[i,j]. Inputs are
;streamline orbit elements and size parameter a,e,w,beta, ar=1, dust inclinations
;I_dust (in radians), and dust size distribution q_size. Set display_orbits=1
;to plot each pair of crossed orbits within distance display_range. Outputs are
;alpha_star, as well as a cumulative count of crossing distances N_crossing versus
;distance r_crossing.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro calc_alpha, inputs, constants, a, e, beta, ar, a_ring, ring_idx, alpha_star, $
    r_crossing, N_crossing, N_collisions, N_disruptions, grazing_orbits=grazing_orbits

;Extract relevant inputs
N_longitudes = inputs.N_longitudes
I_dust = inputs.I_dust
q_size = inputs.q_size
c_disk = inputs.c_disk
display_orbits = inputs.display_orbits
display_range = inputs.display_range
delay = inputs.delay
Q_star_cgs = inputs.Q_star_cgs
M_star_solar = inputs.M_star_solar
ar_au = inputs.ar_au

;Extract constants
G_cgs = constants.G_cgs
M_sun_gm = constants.M_sun_gm
cm_per_au = constants.cm_per_au

;Number of dust size bins
N_sizes = n_elements(a)

;Max beta
beta_max = max(beta)

;Dust periapse distances
q_peri = a*(1 - e)

;Storage for alpha_star array
alpha_star = fltarr(N_sizes, N_sizes)

;Orbital speed^2 of outer planetesimal ring, in (cm/sec)^2
M_star_gm = M_star_solar*M_sun_gm
ar_cm = ar_au*cm_per_au
V2_out_cgs = G_cgs*M_star_gm/ar_cm

;Storage for cumulative count of crossing distances
Nc = 1001
r_min = 0.99*min(q_peri)
r_max = 5.0
r_crossing = (r_max - r_min)*findgen(Nc)/(Nc - 1) + r_min
N_crossing = fltarr(Nc)

;Storage for number of collisions and disruptions experienced by each streamline
N_collisions = lonarr(N_sizes)
N_disruptions = lonarr(N_sizes)

;Impactor longitude of periapse *relative* to target's longitude of periapse.
twopi = 2.0*!pi
delta_w = twopi*findgen(N_longitudes)/N_longitudes

;For plotting orbits
N_plot = 10001
theta = !pi*(2*findgen(N_plot)/(N_plot - 1) - 1)
rng = display_range*[ -1, 1 ]
print, format = '($, A0)', 'Calculating alpha_star: '

;The default grazing_orbits=0 is set so that the code will ignore grazing orbits.
;Grazing orbits are orbits that run tangent at one site rather than
;crossing at two sites. Grazing orbits have large
;collision probabilities, and allowing grazing orbits generates dust abundances
;that are `ragged' when plotted versus dust size. Ignoring grazing orbits
;gets rid of the raggedness, but generates results that are otherwise similar.
;Set grazing_orbits = 1 to account for grazing orbits when
if (keyword_set(grazing_orbits) eq 0) then  grazing_orbits = 0

;Loop over all target dust sizes
for i = 0, N_sizes - 1 do begin $

    ;Loop over all impactor dust sizes
    for j = 0, N_sizes -1 do begin $

        ;Loop over all possible *relative* longitudes of periapse delta_w
        for kj = 0, N_longitudes - 1 do begin $

            ;orbit elements and grain size for target orbit i
            ai = a[i] & $
            ei = e[i] & $
            wi = 0.0 & $             ;this is arbitrary, so set it to zero
            beta_i = beta[i] & $
            ring_idx_i = ring_idx[i] & $

            ;orbit elements and grain size for impactor orbit j
            aj = a[j] & $
            ej = e[j] & $
            wj = delta_w[kj] & $    ;impactor's longitude of peri relative to target's
            beta_j = beta[j] & $
            ring_j = a_ring[j] & $
            ring_idx_j = ring_idx[j] & $

            ;Note that dust grains emitted by the same planetesimal ring with the same
            ;longitudes of periapse (ie, when kj=0) have zero relative velocity
            ;where they intersect, and don't shatter. So proceed only if grains
            ;i and j have a chance of colliding and shattering.
            if ( (ring_idx_i ne ring_idx_j) or (kj ne 0) ) then begin $

                ;find polar coordinates (rc, theta_c) where orbits cross
                orbit_crossing, ai, ei, wi, aj, ej, wj, r_c, theta_c & $

                ;Count the number of sites N_cross where these two orbits cross.
                ;Orbits don't cross if N_cross=0, are grazing if N_cross=1, and
                ;cross at two sites when N_cross=2. Proceed depending on whether
                ;grazing orbits are or are not allowed.
                idx = where(r_c gt 0, N_cross) & $
                proceed = 0 & $
                if ((grazing_orbits eq 0) and (N_cross ge 2)) then proceed = 1 & $
                if ((grazing_orbits eq 1) and (N_cross ge 1)) then proceed = 1 & $
                if (proceed eq 1) then begin $

                    ;Loop over both crossing sites, and get velocity of target (index i)
                    ;& impactor (index j) at collision site r_c[l] in units of ring's
                    ;speed.
                    for l = 0, 1 do begin $
                        if (r_c[l] gt 0.0) then begin $

                            ;include r_c in histogram of the crossing sites, & count 
                            ;collision
                            idx = where(r_crossing ge r_c[l], N_idx) & $
                            if (N_idx gt 0) then $
                                N_crossing[idx[0]] = N_crossing[idx[0]]+1 & $
                            N_collisions[i] = N_collisions[i] + 1 & $

                            ;Get velocity of target i and impactor j at crossing site,
                            ;in units of the outer planetesimal ring's orbital speed.
                            ;Then calculate grains' relative speed^2, in (cm/sec)^2,
                            ;and sin_phi
                            velocity, ai, ei, wi, beta_i, r_c[l], theta_c[l], ar, $
                                I_dust, v_r_i, v_theta_i, v_z_i, v_i & $
                            velocity, aj, ej, wj, beta_j, r_c[l], theta_c[l], ar, $
                                I_dust, v_r_j, v_theta_j, v_z_j, v_j & $
                            v2_rel = (v_r_i - v_r_j)^2 + (v_theta_i - v_theta_j)^2 + $
                                v_z_i^2 + v_z_j^2 & $
                            v2_rel_cgs = v2_rel*V2_out_cgs & $
                            vi_cross_vj = abs( v_r_i*v_theta_j - v_theta_i*v_r_j ) & $
                            sin_phi = vi_cross_vj/v_i/v_j & $

                            ;Print warning if sin_phi=0, which would be problematic 
                            ;because alpha_star then diverges. The happened rarely
                            ;while this code was being developed, and not at all
                            ;for the simulations described in the 2010 ApJ paper.
                            if (sin_phi eq 0) then $
                                print, string(10B) + 'CALC_ALPHA(): WARNING    i = ' $
                                    + string(i) + '    j = ' + string(j) + '    wj = ' $
                                    + string(wj) + '    sin_phi = ' + string(sin_phi) & $

                            ;Proceed if sin_phi>0
                            if (sin_phi gt 0) then begin $

                                ;Grains i and j fragment if disrupt > 0.
                                disrupt = 0.5*v2_rel_cgs/Q_star_cgs - $
                                    ((beta_i^3 + beta_j^3)^2)/((beta_i*beta_j)^3) & $
                                if (disrupt gt 0) then begin & $

                                   ;Add contribution to alpha_star[i,j]
                                    reduced_beta = 1/beta_i + 1/beta_j & $
                                    denominator = N_longitudes*4*sin_phi $
                                        *sqrt(2*aj/r_c[l]-1) & $
                                    alpha = (beta_max^2)*sqrt(1 - beta_i) $
                                        *(reduced_beta^2)*((beta_j/beta_max)^q_size) $
                                        *(ar/r_c[l])*(ar/aj)*((ar/ai)^1.5) $
                                        *((ring_j/ar)^c_disk)/denominator & $
                                    alpha_star[i, j] = alpha_star[i, j] + alpha & $

                                    ;Count fragmentation
                                    N_disruptions[i] = N_disruptions[i] + 1 & $

                                endif & $

                            endif & $

                         endif & $
                     endfor & $

                     ;plot orbits, if desired
                     if (display_orbits eq 1) then begin $
                         vi = theta - wi & $
                         ri = ai*(1.0 - ei^2)/(1.0 + ei*cos(vi)) & $
                         vj = theta - wj & $
                         rj = aj*(1.0 - ej^2)/(1.0 + ej*cos(vj)) & $
                         plot, ri, theta, xrange=rng, yrange=rng, polar=1, xstyle=1, $
                             ystyle=1, title='i = ' + strtrim(string(i), 2) + $
                             '    j = ' + strtrim(string(j), 2) + '    kj = ' + $
                             strtrim(string(kj), 2), charsize=1.5 & $
                         oplot, rj, theta, polar=1, color=128 & $
                         plots, [0, ai*(1.0 - ei)], [0, 0] & $
                         mn = min(rj, idx) & $
                         plots, [0, rj[idx]*cos(theta[idx])], $
                             [0, rj[idx]*sin(theta[idx])], color=128 & $
                         plots, 0.0, 0.0, psym=1 & $
                         l = where(r_c gt 0.0, Nl) & $
                         if (l[0] ne -1) then $
                             oplot, r_c[l], theta_c[l], psym=8, polar=1, symsize=1.5 & $
                         wait, delay & $
                     endif & $

                 endif & $

             endif & $

         endfor & $
     endfor & $
     if ( (i mod ((N_sizes/10)>1)) eq 0) then begin $
         fraction = (100l*i)/N_sizes & $
         print, format = '($, A0)', strtrim(string(fraction),2)+'% ' & $
     endif & $
endfor
print, '100%'

;Convert histogram into cumulative count of N_crossing versus r_crossing
for j = 1, Nc - 1 do N_crossing[j] = N_crossing[j-1] + N_crossing[j]
if (N_crossing[Nc-1] gt 0) then N_crossing = N_crossing/N_crossing[Nc-1]

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Runga Kutta integrator solves for dust populations n versus time tau,
;where n[i, tm] = population  of streamline i at time tau[tm]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro evolve_system, alpha, N_times, N_longitudes, tau1, tau2, tau_stop, n, tau

;Dimensionless time axis tau is logarithmic over times tau1 < tau < tau2, with tau[0]=0
tau = ( alog(tau2) - alog(tau1) )*findgen(N_times - 1)/(N_times - 2) + alog(tau1)
tau = [ 0.0, exp(tau) ]

;Number of streamlines
sz = size(alpha)
N_sl = sz[1]

;Initial n(tau=0) = 0
n0 = fltarr(N_sl)

;Try this timestep first
dtau_try = (1.0e-2)*tau[1]

;Warn when RK integrator chooses a timestep smaller than this threshold
dtau_min = (1.0e-8)*tau[1]

;Name of procedure that calculates derivative dn/dtau
pro_name = 'dd_derivs'

;Fractional error allowed in the Runga Kutta integration of this system over time
rk_error  = 1.0e-5

;Create a structure used to pass alpha and tau_stop parameters to dd_derivs().
params = {alpha:alpha, tau_stop:tau_stop}

;Perform Runge Kutte integration to solve this system of differential equations.
print, format = '($, A0)', 'Evolving system over time:     '
int_sys_eqns, tau, n0, rk_error, dtau_try, dtau_min, pro_name, n, params=params

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Calculate the derivatives dn_dt[j] = derivative dn/dt at time t=tau for streamline j.
;This proceedure is called by subroutines in the int_sys_eqns proceedure.
;Structure params passes the alpha array and tau_stop to this procedure. This procedure
;uses matrix math to avoid any explicit summations, which should make it fast (I think).
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro dd_derivs, t, n, dn_dt, params=params

;extract input parameters
alpha = params.alpha
tau_stop = params.tau_stop

;Calculate dn/dt[beta] due to collisions only
dn_dt = (-n)*(alpha # n)

;Add dust production term if t<tau_stop or tau_stop<0
if ((t lt tau_stop) or (tau_stop lt 0)) then begin $

    dn_dt = dn_dt + 1 & $

    ;Zero elements in |dn_dt| in dn/dt that are sufficiently small,
    ;which speeds up execution times considerably
    dn_dt_max = 1.0e-5 & $
    j = where(abs(dn_dt) lt dn_dt_max) & $
    if (j[0] ne -1) then dn_dt[j] = 0.0 & $

endif

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;This function is called by mpfitfun, which is called by the equilibrium proceedure;
;this function returns the right-hand side of Eqn (27) in Hahn (2010). IDL's matrix math
;is used here, which purportedly speeds up the calculation. The alpha array is
;passed using the _EXTRA parameter that is described in mpfitfun.pro.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function equilibrium_fn, x, n, _EXTRA=extra

alpha = extra.alpha
eqlbrm = n*(alpha # n) - 1

return, eqlbrm
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Calculate dust equilibrium abundances n[beta], which is the solution to Eqn (27)
;in Hahn (2010) at time t^star->infinity. Because dn/dt=0 when in equilibrium, the
;differentials dissapear, which results in a system of nonlinear equations that
;are solved using Craig Markwardt's mpfitfun algorithm.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro equilibrium, alpha, n

print, 'Calculating equilibrium dust abundances...'
sz = size(alpha)
N_sl = sz[1]

;Some dummy variables required by mpfitfun.
x = findgen(N_sl)
y = fltarr(N_sl)

;Allowed tolerance in the solution
err = fltarr(N_sl) + 1.0e-3

;Initial guess is just 1
n_start = fltarr(N_sl) + 1.0

;Limit values of possible solutions such that 0<n<100. 
parinfo = replicate({limited:[1, 1], limits:[0.0, 100.0]}, N_sl)

;Solve for n[beta] that satisfies Eqn (27)
n = mpfitfun('equilibrium_fn', x, y, err, n_start, params=params, quiet=1, $
    functargs = {alpha:alpha}) 

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Calculate dust optical depth tau versus distance r in units of outermost planetesima
;ring radius ar. Inputs are number of grains N in each beta bin, orbit elements a,e, size beta,
;smin_ar=smallest grain radius in units of ring radius, and ring radius ar=1.
;Output is optical depth tau[r] and A_enc=total dust cross section interior to r,
;in units of ring area.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro optical_depth, N, a, e, beta, smin_ar, ar, r, tau, A_enc

;number of different grain sizes
N_sizes = n_elements(a)

;Maximum beta
beta_max = max(beta)

;Storage for enclosed area A_enc[r]
N_r = n_elements(r)
A_enc = dblarr(N_r)
twopi = 2.0*!pi

;loop over all streamlines, and add their contribution to A_enc[r]
for i = 0l, N_sizes - 1l do begin $
    cos_Ecc = (1.0 - r/a[i])/e[i] & $          ;cosine(Eccentric anomaly)
    cos_Ecc = cos_ecc>(-1.0)<1.0 & $           ;this zeros contribution from streamlines
    Ecc = acos(cos_Ecc) & $                    ;               that don't cross radius r
    t_over_T = (Ecc - e[i]*sin(Ecc))/twopi & $ ;time since peri-passage/orbit period
    A_enc = A_enc + 2.0*(N[i]*smin_ar)*smin_ar*t_over_T*((beta_max/beta[i])^2) & $
endfor

;Calculate optical depth
dA = deriv(r, A_enc)
tau = 0.5*dA/r

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Calculate the disk's total cross section and mass. Inputs are the number of dust
;grains N[beta] in each beta size-bin, smin_ar=grain radius in units of ring radius,
;ar_km=ring radius in km, mass_min_gm=mass of smallest grain in grams. Outputs are
;A_disk_ringarea=total dust cross section in units of the ring's area, A_disk_km2
;(units of km^2), M_disk_smallestgrain=total dust mass in units of smallest grain
;mass, and M_disk_gm=dust mass in gm.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro disk_totals, N, beta, smin_ar, ar_km, mass_min_gm, A_disk_ringarea, A_disk_km2, $
    M_disk_smallestgrain, M_disk_gm

;Number of distinct dust sizes and max beta
N_sizes = n_elements(beta)
beta_max = max(beta)

;Total disk cross section in units of the ring area Pi*ar^2, and km^2
A_disk_ringarea = total( (N*smin_ar)*smin_ar*((beta_max/beta)^2) )
A_disk_km2 = A_disk_ringarea*!pi*(ar_km^2)

;Total disk mass in units of smallest grain's mass, and grams
M_disk_smallestgrain = total( N*((beta_max/beta)^3) )
M_disk_gm = M_disk_smallestgrain*mass_min_gm

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Henyey-Greenstein phase function. Inputs are the light-scattering-asymmetry-factor g,
;and the scattering angle scat_angle in radians. This function returns the value of
;the phase function in units of 1/(one steradian).
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function hg_phase_function, g, scat_angle

factor = (1 - g^2)/(4*!pi)
den = 1 + g^2 - 2*g*cos(scat_angle)
phs_fn = factor/(den^1.5)

return, phs_fn
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Calculate edge-on disk's total surface brightness B (in unfiltered light) versus
;projected distance ell. B has units of B0, where that constant has units of
;ergs/cm^2/sec/steradian (B0_ster) and ergs/cm^2/sec/arcsec^2 (B0_arc2), while ell
;has units of disk radius ar. To get the disk's surface brightness in a narrowband
;filter, use B*B0_arc2_cm, which has units of ergs/cm^2/sec/arcsec^2/cm.
;Also, setting N_longitudes_alternate>N_longitudes will give a smoother result at
;large x.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro surface_brightness, a, e, beta, N, x, ar, inputs, constants, smin_ar, B, $
    B0_cgs_per_ster, B0_cgs_per_arc2, B0_Jansky_per_arc2, lambda_wien_microns, $
    N_longitudes_alternate=N_longitudes_alternate

print, 'Calculating disk surface brightness...'

;Extract inputs
N_longitudes = inputs.N_longitudes
I_dust = inputs.I_dust
g_dust = inputs.g_dust
Q_scat = inputs.Q_scat
ar_au = inputs.ar_au
wavelength_nm = inputs.wavelength_nm
T_star_K = inputs.T_star_K
L_star_solar = inputs.L_star_solar
phase_law = inputs.phase_law

;Extract constants
cm_per_au = constants.cm_per_au                ;number of cm in one AU
nm_per_cm = constants.nm_per_cm                ;number of nm in one cm 
micron_per_cm = constants.micron_per_cm        ;number of microns in one cm 
arc2_per_ster = constants.arc2_per_ster        ;number of arcseconds in one steradian
L_sun_cgs = constants.L_sun_cgs                ;solar luminosity in cgs unigs
c_cgs = constants.c_cgs                        ;speed of light in cgs units
sbc_cgs = constants.sbc_cgs                    ;stefan-boltzmann constant in cgs
Jansky = constants.Jansky                      ;one Jansky in cgs units

;Adjust N_longitudes as needed
if (keyword_set(N_longitudes_alternate) eq 1) then N_longitudes = N_longitudes_alternate

;Number of distinct streamlines
N_sizes = n_elements(a)

;Number of elements in ell axis
N_x = n_elements(x)

;Max beta
beta_max = max(beta)

;Star's luminosity in cgs units
L_star_cgs = L_star_solar*L_sun_cgs

;Ring semimajor axis in cm
ar_cm = ar_au*cm_per_au

;Constant B0 in ergs/cm^2/sec/steradian (B0_cgs_per_ster),
;and ergs/cm^2/sec/arcsec^2 (B0_cgs_per_arc2)
Omega_1 = 1.0                   ;solid angle of 1 steradian, in units of steradians
B0_cgs_per_ster = Q_scat*L_star_cgs/(16d*!pi*I_dust*Omega_1)/ar_cm/ar_cm
B0_cgs_per_arc2 = B0_cgs_per_ster/arc2_per_ster

;Frequency (in Hz) of narrowband observation from observation wavelength (nanometers).
wavelength_cm = wavelength_nm/nm_per_cm                        ;wavelength in cm
frequency_Hz = c_cgs/wavelength_cm

;Evaluate Pi*Planck function of frequency, in cgs units (ergs/cm^2/sec/Hz)
planck_fn, T_star_K, frequency_Hz, constants, planck_fn_cgs

;For narrowband observations, multiply B0 by factor f_n = Pi*Planck function/(SBC*T^4)
;(which has units of 1/Hz) to get surface brightness narrowband observations
f_n = planck_fn_cgs/sbc_cgs/(T_star_K^4)

;Constant B0 in units of Janskies/arcsec^2
B0_Jansky_per_arc2 = B0_cgs_per_arc2*f_n/Jansky

;Wavelength of star's peak emission from Wien's law
lambda_wien_cm = 0.28978/T_star_K
lambda_wien_microns = micron_per_cm*lambda_wien_cm

;Storage for disk's edge-on surface brightness, in units of constant B0
B = dblarr(N_x)

;Loop over all streamlines, find the two sites where the line-of-sight (LOS) intersects
;the orbit, and add that streamline's contribution to disk's surface brightness B[x]
for i = 0l, N_sizes - 1l do begin $            ;loop over all sizes beta[i]   
    for j = 0l, N_longitudes - 1l do begin $   ;loop over longitudes of periapse w[j]

        ;The streamline's longitude of periapse
        w = (2.0*!pi*j)/N_longitudes & $

        ;calculate x_max = max(x) where LOS will intersects orbit i,j
        esinw = e[i]*sin(w) & $
        ecosw = e[i]*cos(w) & $
        p = a[i]*(1 - e[i]^2) & $
        factor = sqrt( 1 - (esinw)^2 ) & $
        x_max = p*factor/( 1 + ecosw*factor - (esinw)^2 ) & $

        ;These LOSs x[k] intersect the ellipse
        k = where(x lt x_max, N_k) & $

        ;calculate r1(x) and r2(x) = two radial distances where LOS intersects orbit i,j
        if (N_k gt 0) then begin $
            Ax = p*cos(w) - e[i]*x & $
            Bx = p*sin(w) & $
            phi = atan(Bx, Ax) & $
            arg = x/sqrt( Ax^2 + Bx^2) & $
            angle = fltarr(N_x) & $
            angle[k] = acos(arg[k]) & $
            v1 =  angle - phi & $
            v2 = -angle - phi & $
            theta1 = v1 + w & $
            theta2 = v2 + w & $
            r1 = p/(1 + e[i]*cos(v1)) & $
            r2 = p/(1 + e[i]*cos(v2)) & $
            r_k = [ [r1], [r2] ] & $
            theta_k = [ [theta1], [theta2] ] & $
        endif & $

        ;Set check=1 to plot orbit, LOS, and intersection sites
        check = 0 & $
        if (check eq 1) then begin $
            N_plot = 1001 & $
            v_plot = 2.0*!pi*findgen(N_plot)/(N_plot - 1) & $
            r_plot = p/(1 + e[i]*cos(v_plot)) & $
            theta_plot = v_plot + w & $
            rng = 1.05*max(r_plot)*[-1, 1] & $
            plot, r_plot, theta_plot, xrange=rng, yrange=rng, xstyle=1, ystyle=1,  $
                title = 'i,j = ' + string(i)+' '+string(j), polar=1 & $
            plots, 0, 0, psym=1 & $
            plots, [x_max, x_max], rng, linestyle=2 & $
            for idx = 0, N_k - 1 do begin $
                if (x[k[idx]] lt x_max) then begin $
                    x_dot = r_k[k[idx], *]*cos(theta_k[k[idx], *]) & $
                    y_dot = r_k[k[idx], *]*sin(theta_k[k[idx], *]) & $
                    plots, x_dot, y_dot, psym=8 & $
                    oplot, x[k[idx]]*[1, 1], rng, color=128 & $
                endif & $
            endfor & $
            strng = '' & $
            read, strng & $
        endif & $

        ;loop over sites where LOS intersects orbit, and add brightness contribution
        if (N_k gt 0) then begin $
            for idx = 0, 1 do begin $
                r = r_k[k, idx] & $
                theta = theta_k[k, idx] & $
                ;evaluate phase function x 1 steradian (which is dimensionless)
                scat_angle = 0.5*!pi + theta & $
                if (phase_law eq 'HG') then $
                    phase_fn = hg_phase_function(g_dust, scat_angle) & $
                den = abs( sin(theta) + e[i]*sin(w) ) & $
                B_k = phase_fn*(N[i]*smin_ar)*smin_ar*sqrt(1 - e[i]^2)*(ar/a[i]) $
                    *((ar/r)^3)*((beta_max/beta[i])^2)/N_longitudes/den & $
                B[k] = B[k] + B_k & $
            endfor & $
        endif & $

    endfor & $
endfor

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Evaluate Pi*Planck blackbody function of *frequency*. Inputs are the effective
;temperature (in units of Kelvin), and frequency (in Hz). Output is Pi*Planck
;function of frequency in cgs units (ergs/cm^2/sec/Hz).
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro planck_fn, T_K, frequency_Hz, constants, planck_fn_cgs

;Extract physical constants
c_cgs = constants.c_cgs                        ;speed of light in cgs units
h_cgs = constants.h_cgs                        ;Planck's constant in cgs units
k_cgs = constants.k_cgs                        ;Boltzmann constant in cgs units
nm_per_cm = constants.nm_per_cm                ;number of nm in one cm 

;Evaluate Pi*Planck function, in cgs units (ergs/cm^2/sec/Hz)
argument = h_cgs*frequency_Hz/(k_cgs*T_K)
Planck_fn_cgs = ( 2.0*!pi*h_cgs*(frequency_Hz^3)/(c_cgs^2) )/( exp(argument) - 1.0 )

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

