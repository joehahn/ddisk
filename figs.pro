;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;figs.pro
;
;Read the output generated by v4.0 of ddisk.pro, do some additional calculations,
;and make some figures
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Set pflag=1 to generate some postscript plots
pflag=0

;Read inputs and outputs, and compile subroutines.pro
.run subroutines
restore, filename = 'ddisk.dat'
@inputs

;Extract some constants
micron_per_cm = constants.micron_per_cm
cm_per_au = constants.cm_per_au
micron_per_cm = constants.micron_per_cm
sec_per_yr = constants.sec_per_yr
km_per_au = constants.km_per_au
M_moon_gm = constants.M_moon_gm
M_earth_gm = constants.M_earth_gm

;Radiation pressure efficiency
Q_rp = 1 - g_dust*Q_scat

;Re-calculate dust radii in microns and masses in gm
beta_max = max(beta)
s_min_microns = 0.58*Q_rp*L_star_solar/M_star_solar/dust_density_cgs/beta_max
s_microns = (beta_max/beta)*s_min_microns
s_cm = s_microns/micron_per_cm
mass_gm = (4.0*!pi/3.0)*dust_density_cgs*(s_cm^3)
mass_min_gm = min(mass_gm)

;ring semimajor axis in cm and microns
ar_cm = ar_au*cm_per_au
ar_microns = ar_cm*micron_per_cm

;Smallest grain radius in units of ring semimajor axis
smin_ar = s_min_microns/ar_microns

;Number of distinct streamlines
N_sizes = n_elements(a)

;Ring orbit period in years and seconds
Tr_yr = sqrt( (ar_au^3)/(M_star_solar) )
Tr_sec = Tr_yr*sec_per_yr

;S_qc and p_max (in grains/sec) are related to the dust mass production rate
exponent = q_size - 3.0
factor = ((beta/beta_max)^exponent)*((a_ring/ar)^c_disk)
S_qc = total(factor)
P0 = double(M_dpr_cgs)/S_qc/mass_min_gm
print, 'dust mass production rate (gm/sec)             = ', M_dpr_cgs
print, 'execution times (minutes)                      = ', execution_time_min
print, 'number of dust sizes                           = ', N_sizes
print, 'number of longitudes                           = ', N_longitudes
print, 'number of planetesimal rings                   = ', N_rings
print, 'S_qc                                           = ', S_qc
print, 'P0 (grains/sec)                                = ', P0

;Dust production rate P_i, in grains per sec
P_i = P0*((beta/beta_max)^q_size)*((a_ring/ar)^c_disk)

;Plot streamlines' apocenter distances versus beta.
plot_io, beta, Q_apo, psym=8, charsize=1.5, xtitle='!7b!3', $
    ytitle='Q!lapo!n    (a!lr!n)'
u = uniq(ring_idx)
for j = 0, n_elements(u) - 1 do begin $
    idx = where(ring_idx eq ring_idx[u[j]]) & $
    oplot, beta[idx], Q_apo[idx] & $
    oplot, beta[idx], Q_peri[idx], color=128 & $
endfor
print, 'range of dust apocenter distances (ring radius)= ', min(Q_apo), max(Q_apo)
print, 'range of dust beta values                      = ', min(beta), max(beta) 
print, 'range of dust radii (microns)                  = ', min(s_microns),max(s_microns)

;;Check distribution of crossing sites
plot, r_crossing, N_crossing, psym=10, xrange=[0.9, 3], xstyle=1, charsize=1.5, $
    xtitle='distance r/a!lr!n', ytitle='Cumulative count of collision radii r/a!lr!n'
idx = where(N_crossing ge 0.5)
r_50 = r_crossing[idx[0]]
plots, r_crossing[idx[0]], N_crossing[idx[0]], psym=8, symsize=1.5
idx = where(N_crossing ge 0.95)
r_95 = r_crossing[idx[0]]
plots, r_crossing[idx[0]], N_crossing[idx[0]], psym=8, symsize=1.5
print, '50% of dust collisions occur w/in (ring radius)=  ', r_50
print, '95% of dust collisions occur w/in (ring radius)=  ', r_95

;Plot disruption fraction
disrupt_frac = float(n_disruptions)/N_collisions
plot, beta, disrupt_frac, charsize=1.5, xtitle='!7b!3', ytitle='disruption fraction', $
    psym=3
N_collisions_total = total(N_collisions)
N_disruptions_total = total(N_disruptions)
disrupt_frac_total = float(N_disruptions_total)/N_collisions_total
print, 'fraction of collisions that disrupt            =  ', disrupt_frac_total

;If evolve_system='NO', then set n_star[i,*] = n_star_eq[i] at all times,
;to prevent hiccups elsewhere in this code
if (evolve_system eq 'NO') then begin $
    tau = ( alog(tau2) - alog(tau1) )*findgen(N_times - 1)/(N_times - 2) + alog(tau1) &$
    tau = [ 0.0, exp(tau) ] & $
    n_star = fltarr(N_sizes, N_times) & $
    for i = 0, N_sizes - 1 do n_star[i, *] = n_star_eq[i] & $
endif

;Plot n_star[beta, tau] for all streamlines, where n_star=streamlines' rel' abundances
j = where(n_star gt 0.0)
n_rng = [ 0.5*min(n_star[j]), 2.0*max(n_star[j])]
j = where(tau gt 0.0)
tau_rng = [ 0.5*min(tau[j]), 4.0*max(tau[j])]
plot_oo, tau, n_star[0, *], yrange=n_rng, xstyle=1, ystyle=1, xrange=tau_rng, $
    nodata=1, charsize=1.5, xtitle='time !7s!3', ytitle='n(!7b,s!3)'
for i = 0, N_sizes - 1 do begin $
    oplot, tau, n_star[i, *] & $
    xyouts, 1.2*tau[N_times - 1], 0.8*n_star[i, N_times - 1], $
        strtrim(string(beta_idx[i]),2) & $
    plots, tau[N_times - 1], n_star_eq[i], psym=8, color=128 & $
endfor

;Timescale T0 in seconds and years
T0_sec = sqrt(I_dust*Tr_sec/P0)/smin_ar
T0_yr = T0_sec/sec_per_yr
T0_yr_check = (56002d0)*sqrt(S_qc/106.5)*sqrt(I_dust/0.1)*sqrt(L_star_solar) $
    *(M_star_solar^(-0.75))*((M_dpr_cgs/1.0e13)^(-0.5))*((ar_au/50.0)^1.75)
print, 'T0 (yrs)                                       = ', T0_yr, T0_yr_check

;time t in years and seconds
t_yr = tau*T0_yr
t_sec = t_yr*sec_per_yr

;Dust abundance factor N0
N0 = double(P0)*T0_sec
N0_check = (2.54d+34)*(Q_rp^(-2.5))*((S_qc/106.5)^(-0.5))*sqrt(I_dust/0.1) $
    *(L_star_solar^(-2.5))*(M_star_solar^2.25)*(dust_density_cgs^2) $
    *sqrt(M_dpr_cgs/1e13)*((ar_au/50.0)^1.75)
print, 'N0                                             = ', N0, N0_check

;Calculate dust abundance N[size, time] and equilibrium N_eq
N = dblarr(N_sizes, N_times)
for tm = 0, N_times - 1 do N[0, tm] = N0*(P_i/P0)*n_star[*, tm]
N_eq = N0*(P_i/P0)*n_star_eq

;N_total[i, j] = total number of same-sized grains beta[i] versus time tau[j], while
;P_total[i] = total dust production rate (in grains/sec) for grains of size beta[i]
s = sort(beta_idx)
beta_idx_s = beta_idx[s]
u = uniq(beta_idx_s)
beta_idx_u = beta_idx_s[u]
Nu = n_elements(beta_idx_u)
N_total = dblarr(Nu, N_times)
beta_total = fltarr(Nu)
a_total = fltarr(Nu)
e_total = fltarr(Nu)
beta_idx_total = intarr(Nu)
s_microns_total = fltarr(Nu)
P_total = dblarr(Nu)
for i = 0, Nu - 1 do begin $
    j = where(beta_idx eq beta_idx_u[i]) & $
    N_total[i, 0] = transpose(total(N[j, *], 1)) & $
    beta_total[i] = beta[j[0]] & $
    a_total[i] = a[j[0]] & $
    e_total[i] = e[j[0]] & $
    beta_idx_total[i] = beta_idx[j[0]] & $
    s_microns_total[i] = s_microns[j[0]] & $
    P_total[i] = total(P_i[j]) & $
endfor

;equilibrium timescale = time where small grains reach two-thirds their max population
mx = max(beta_idx_total, i)
N_max = N_total[i, N_times - 1]                                       
tm = where(N_total[i, *] gt N_max*2.0/3.0)                                        
if (evolve_system eq 'YES') then begin $
    T_eq_yr = t_yr[tm[0]] & $
endif else T_eq_yr = 130.0*T0_yr
N_approx = N_max*tanh(t_yr/T_eq_yr)
print, 'T_eq (Myrs)                                    = ', T_eq_yr/1.0e6
print, 'T_eq/T0                                        = ', T_eq_yr/T0_yr

;Set N_total_final[beta]=dust abundance at stellar age
tm_age = where(t_yr ge t_star_yr)
tm_age = tm_age[0]
if (tm_age[0] eq -1) then tm_age = N_times - 1
N_total_final = N_total[*, tm_age]
N_final = N[*, tm_age]

;Plot N_total versus time for selected streamlines
pflag = 0
tau_rng = [1.0e-4, 1.0e4]
t_rng = tau_rng*T0_yr
NN_rng = [5.0e-7, 5.0e2]
N_rng = NN_rng*N0
betas = [ '0.496', '0.48', '0.45', '0.38', '0.30', '0.22', '0.15', '0.10', '0.07']
thck = 3*pflag + 1
if (pflag eq 1) then setplot
plot, tau, N_total[0, *], yrange=NN_rng, xstyle=9, ystyle=9, xrange=tau_rng, $
    nodata=1, charsize=1.5, xtitle='time t/T!l0!n', ythick=thck, charthick=thck, $
    ytitle='relative dust abundance N!li!n/N!l0!n', thick=thck, xthick=thck, $
    xlog=1, ylog=1
for j = 0, n_elements(betas) - 1 do begin $
    i = where(beta_total ge float(betas[j])) & $
    i = i[0] & $
    if (i ge 0) then begin $
        oplot, tau, N_total[i, *]/N0, thick=thck & $
        xyouts, 120, 1.3*N_total[i, N_times - 1]/N0, '!7b!3 = ' + $
            strmid(betas[j], 0, 6), charsize=1.1, charthick=2*pflag + 1 & $
        plots, tau[tm_age], N_total_final[i]/N0, psym=8 & $
    endif & $
endfor
oplot, tau, N_approx/N0, color=128, thick=thck+6, linestyle=2
;oplot, tau[tm_age] + [0, 0], NN_rng, color=128
axis, yaxis=1, charsize=1.5, ystyle=1, yrange=N_rng, $
    ytitle='!Ctotal dust abundance N!li!n(t)', ythick=thck, charthick=thck
axis, xaxis=1, charsize=1.5, xstyle=1, xrange=t_rng, $
    xtitle='time t    (yrs)!C', xthick=thck, charthick=thck
if (pflag eq 1) then output_plot, 'N.ps'
print, 'tau at age of star                             = ', tau[tm_age]
pflag = 0

;Plot N_final and P versus beta
N_rng = [2.0d28, 2.0d36]
if (pflag eq 1) then setplot
thck = 3*pflag + 1
plot_oo, s_microns_total, N_total_final, xrange=[0.9, 11.0], xstyle=1, charsize=1.5, $
    ystyle=1, yrange=N_rng, xtitle='grain radius R    (!7l!3m)', $
    ytitle='equilibrium abundance N(R)', $
    thick=thck, xthick=thck, ythick=thck, charthick=thck
j = where(s_microns_total lt 6.0)
j = j[0]
power_law = N_total_final[j]*(s_microns_total/s_microns_total[j])^(-6.0)
oplot, s_microns_total, power_law, linestyle=2, thick=thck
j = 0
power_law = N_total_final[j]*(s_microns_total/s_microns_total[j])^(-q_size)
oplot, s_microns_total, power_law, linestyle=1, thick=thck
xyouts, 1.25, 7.0d34, 'N(R)', charsize=1.5, charthick=thck
xyouts, 1.11, 3.0d33, 'R!u-6', charsize=1.5, charthick=thck
xyouts, 1.11, 2.0d31, 'R!u-3.5', charsize=1.5, charthick=thck
if (pflag eq 1) then output_plot,'N_final.ps'

;Collisional lifetime, in seconds and years
t_col_sec = N_total_final/(P_total)
t_col_yr = t_col_sec/sec_per_yr
if (pflag eq 1) then setplot
thck = 3*pflag + 1
plot_io, s_microns_total, t_col_yr, xrange=[0, 10], xstyle=1, charsize=1.5, ystyle=1, $
    yrange=[ 1.0e2, 1.0e8], ytitle='collisional lifetime T!lc!n    (yrs)',  thick=thck,$
    xthick=thck, ythick=thck, xtitle='grain radius R    (!7l!3m)', charthick=thck
j = where(s_microns_total le 4.0)
j = j[0]
t_col_yr_approx = t_col_yr[j]*((s_microns_total/s_microns_total[j])^(-2.4))
j = where(s_microns_total ge 2.0*min(s_microns_total))
oplot, s_microns_total[j], t_col_yr_approx[j], color=128, linestyle=2

;Compare collision timescales to PR-drag timescales
T_pr_sec = 0.5d*(ar_cm^2)*constants.c_cgs $
    /((constants.G_cgs*constants.M_sun_gm)/M_star_solar)
T_pr_yr = T_pr_sec/constants.sec_per_yr
T_a_yr = (T_pr_yr/beta_total)*((a_total/ar)^2)*((1 - e_total^2)^1.5)/(1 + 1.5*e_total^2)
T_e_yr = 0.8*(T_pr_yr/beta_total)*sqrt(1 - e_total^2)*((a_total/ar)^2)
oplot, s_microns_total, T_a_yr, thick=thck
oplot, s_microns_total, T_e_yr, linestyle=2, thick=thck
axis, yaxis=1, charsize=1.5, ystyle=1, ythick=thck, charthick=thck, $
    ytitle='PR drag timescales T!la!n and T!le!n    (yrs)'
xyouts, 1.5, 1.0e5, 'T!lc!n', charsize=1.5, charthick=thck
xyouts, 7.0, 5.0e7, 'T!la!n', charsize=1.5, charthick=thck
xyouts, 7.0, 1.0e7, 'T!le!n', charsize=1.5, charthick=thck
if (pflag eq 1) then output_plot,'T_c.ps'

;Generate radial axis r (in units of ar) over which the optical depth is calculated.
;Keep in mind that tau varies as the reciprocal of the dust radial velocity,
;which is zero at periapse and apoapse, so choose r-values appropriately to avoid
;spurious peaks in opt_depth[r].
r0 = [Q_peri, Q_apo]
r0 = r0[sort(r0)]
r = 0.5*(r0 + shift(r0, -1))
rm1 = abs(r - 1)
j = where(rm1 gt 0.05, Nj)
r = r[j[0:Nj-2]]
s = sort(r)
r = r[s]
u = uniq(r)
r = r[u]
r = [1.02, 1.07, 1.15, r, 40.0, 60.0, 1.01*max(Q_apo) ]
s = sort(r)
r = r[s]

;Normal optical depth opt_depth. Note this curve gets ragged in the disk's outer parts,
;at r~Q_apo_max/2, due to the sparse number of dust orbits that range out to these distances 
optical_depth, N_final, a, e, beta, smin_ar, ar, r, opt_depth, A_enc
x_rng = [0.4, 120]
tau_rng = [ 1.0e-6, 1.0e-2] 
plot_oo, r, opt_depth, xtitle='distance r    (a!lr!n)', charsize=1.5, ystyle=1, $
    ytitle='optical depth !7s!3(r)', xrange=x_rng, xstyle=1, yrange = tau_rng
oplot, r, opt_depth, psym=8
j = where(r gt 8.0)
j = j[0]
oplot, r, opt_depth[j]*(r[j]/r)^(1.5), color=128

;Normal optical depth at selected times
taus = [ 0.01, 0.1,  1.0, 99.9, 100.0, 150.0, 1000.0 ]
x_rng = [0.5, 120]
tau_rng = [ 3.0e-8, 3.0e-3] 
N_taus = n_elements(taus)
if (pflag eq 1) then setplot
thck = 3*pflag + 1
for j = 0, N_taus - 1 do begin $
    tm = where(tau ge taus[j]) & $
    tm = tm[0] & $
    optical_depth, N[*, tm], a, e, beta, smin_ar, ar, r, opt_depth_tau, A_enc_tau & $
    if (j eq 0) then plot_oo, r, opt_depth, xtitle='distance r    (r!lout!n)', $
        charsize=1.5, ystyle=1, ytitle='optical depth !7s!3(r)', xrange=x_rng, $
        xstyle=1, yrange = tau_rng, nodata=1, thick=thck, xthick=thck, ythick=thck, $
        charthick=thck & $
    oplot, r, opt_depth_tau, thick=thck, psym=0 & $
endfor
xp = [9.0, 19.5]
power_law = (0.8e-3)/(xp^1.5)
oplot, xp, power_law, thick=thck
xyouts, 13.5, 1.85e-5, 'r!u-3/2!n', charsize=1.3, charthick=thck
xyouts, 2.77, 2.0e-7, 't!u*!n=0.01', charsize=1.3, charthick=thck
xyouts, 6.2, 4.0e-7, '0.1', charsize=1.3, charthick=thck
xyouts, 11.0, 8.6e-7, '1', charsize=1.3, charthick=thck
xyouts, 14.5, 4.4e-6, '100', charsize=1.3, charthick=thck
if (pflag eq 1) then output_plot,'tau_time.ps'

;Make sure radial optical depth << 1
pflag = 0
integrnd = opt_depth/(2*I_dust*r)
Nr = n_elements(r)
opt_depth_r = fltarr(Nr)
for j = 1, Nr - 1 do opt_depth_r[j] = int_tabulated(r[0:j], integrnd[0:j])
if (pflag eq 1) then setplot
thck = 3*pflag + 1
plot_oo, r, opt_depth_r, xtitle='distance r    (r!lout!n)', charsize=1.5, ystyle=1, $
    ytitle='radial optical depth !7s!3!lr!n(r)', xrange=[0.5, 150], xstyle=1, $
    yrange = [1.0e-4, 1.0e-2], thick=thck, xthick=thck, ythick=thck, charthick=thck & $
if (pflag eq 1) then output_plot,'radial_opt_depth.ps'
pflag = 0

;Verify that enclosed dust cross section A_enc[r] = integral(2*r*optical depth)
integrand = 2.0*r*opt_depth
N_r = n_elements(r)
A_enc_check = fltarr(N_r)
for j = 1, N_r - 1 do A_enc_check[j] = int_tabulated(r[0:j], integrand[0:j])
plot_oo, r, A_enc, xrange=[1,210], xtitle='distance r    (a!lr!n)', charsize=1.5, $
    ytitle='enclosed cross section A!lenc!n(r)/!7p!3a!lr!n!u2!n', xstyle=1
oplot, r, A_enc_check, color=128

;total dust cross section and mass
ar_km = ar_au*km_per_au
disk_totals, N_final, beta, smin_ar, ar_km, mass_min_gm, A_disk_ringarea, A_disk_km2, $
    M_disk_smallestgrain, M_disk_gm
M_disk_moon = M_disk_gm/M_moon_gm
M_disk_earth = M_disk_gm/M_earth_gm
M_disk_lost_earth = M_dpr_cgs*t_star_yr*sec_per_yr/M_earth_gm
print, 'total dust cross section (ring area)           = ', A_disk_ringarea, max(A_enc)
print, 'total dust cross section (km^2)                = ', A_disk_km2
print, 'total dust mass (gm)                           = ', M_disk_gm
print, 'total dust mass (lunar masses)                 = ', M_disk_moon
print, 'total dust mass (earth masses)                 = ', M_disk_earth
print, 'age of star (Myrs)                             = ', t_star_yr/1.0e6
print, 'total mass lost over star lifetime (earthmass) = ', M_disk_lost_earth

;Calculate edge-on disk's surface brightness B[ell], which is in units of 
;egs/cm^2/sec/steradian (B0_cgs_per_ster) or egs/cm^2/sec/arcsecond^2 (B0_cgs_per_arc2)
;in unfiltered light, or in Janskies per arcsec^2 (B0_Jansky_per_arc2) for a 
;narrowband observation. Note that this curve can seem ragged when the line-of-site
;passes near an orbit's periapse or apoapse.
Nx = 101
x_min = 0.15
x_max = 200
x = (alog(x_max) - alog(x_min))*findgen(Nx)/(Nx - 1) + alog(x_min)
x = exp(x)
surface_brightness, a, e, beta, N_final, x, ar, inputs, constants, smin_ar, B, $
    B0_cgs_per_ster, B0_cgs_per_arc2, B0_Jansky_per_arc2, lambda_wien_microns, $
    N_longitudes_alternate=300
sb_cgs_per_arc2 = B*B0_cgs_per_arc2
plot_oo, x, sb_cgs_per_arc2, xrange=[0.1, 50], xstyle=1, yrange=[2.0e-19, 2.0e-12], $
    ystyle=1, xtitle='projected distance from star x/a!lr!n', charsize=1.5, psym=0,$
    ytitle='disk surface brightness B(x)    (ergs/cm!u2!n/sec/arcsec!u2!n)'
j = where(x gt 8.0)
j = j[0]
oplot, x, sb_cgs_per_arc2[j]*(x[j]/x)^(3.5), color=128
print, 'radius of smallest dust grain (microns)         =  ', min(s_microns)
print, 'wavelength of star`s peak emission (microns)    =  ', lambda_wien_microns

;Surface brightness at selected times
taus = [ 0.01, 0.1,  1.0, 100.0 ]
x_rng = [0.8, 50]
sb_rng = [1.0e-21, 1.0e-12]
N_taus = n_elements(taus)
if (pflag eq 1) then setplot
thck = 3*pflag + 1
for j = 0, N_taus - 1 do begin $
    tm = where(tau ge taus[j]) & $
    tm = tm[0] & $
    surface_brightness, a, e, beta, N[*, tm], x, ar, inputs, constants, smin_ar, B_tau, $
        B0_cgs_per_ster, B0_cgs_per_arc2, B0_Jansky_per_arc2, lambda_wien_microns, $
    N_longitudes_alternate=300 & $
    sb = B_tau*B0_cgs_per_arc2 & $
    if (j eq 0) then plot_oo, x, sb, xrange=x_rng, xstyle=1, charsize=1.5, nodata=1, $
        yrange=sb_rng, ystyle=1, xtitle='projected distance from star x    (r!lout!n)', $
        ytitle='surface brightness B(x)    (ergs/cm!u2!n/sec/arcsec!u2!n)', $
        thick=thck, xthick=thck, ythick=thck, charthick=thck & $ & $
    oplot, x, sb, thick=thck & $
endfor
xp = [6.0, 9.5]
power_law = (6.5e-13)/(xp^3.5)
oplot, xp, power_law, color=128, thick=thck
xyouts, 7.8, 6.4e-16, 'x!u-3.5!n', color=128, charsize=1.3, charthick=thck
xp = [2.2, 3.3]
power_law = 0.25*(1.5e-13)/(xp^5.0)
oplot, xp, power_law, color=128, thick=thck
xyouts, 2.1, 1.5e-16, 'x!u-5.0!n', color=128, charsize=1.3, charthick=thck
xyouts, 2.83, 3.2e-17, 't/T!l0!n=0.01', charsize=1.3, charthick=thck
xyouts, 6.4, 3.2e-17, '0.1', charsize=1.3, charthick=thck
xyouts, 9.3, 3.2e-17, '1', charsize=1.3, charthick=thck
xyouts, 9.3, 7.3e-17, '100', charsize=1.3, charthick=thck
if (pflag eq 1) then output_plot,'sb_time.ps'

;plot filtered surface brightness
sb_Jansky_per_arc2 = B*B0_Jansky_per_arc2
plot_oo, x, sb_Jansky_per_arc2, xrange=[0.5, 50], xstyle=1, yrange=[1.0e-10, 2.0e-4], $
    ystyle=1, xtitle='projected distance from star x/a!lr!n', charsize=1.5, $
    ytitle='disk surface brightness B(x)    (Jy/arcsec!u2!n)'
j = where(x gt 8.0)
j = j[0]
oplot, x, sb_Jansky_per_arc2[j]*(x[j]/x)^(3.5), color=128

;Plot the fractional errors in numerical solution for n_star and N.
;Set check_errors=0 to skip this step.
check_errors = 1
max_residual_n_star = -1.0
max_residual_N = -1.0
if (check_errors eq 1) then begin $
    if (evolve_system eq 'YES') then begin $

        ;Calculate residual in equation for scaled dust abundance n_star.
        print, 'Checking residual in n_star...' & $
        dn_dtau = fltarr(N_sizes, N_times) & $
        for i = 0, N_sizes - 1 do dn_dtau[i, 0] = deriv(tau, n_star[i, *]) & $
        residual_n_star = fltarr(N_sizes, N_times) & $
        for tm = 0, N_times - 1 do begin $
            for i = 0, N_sizes - 1 do begin $
                residual_n_star[i, tm] = dn_dtau[i, tm] - 1 + n_star[i, tm]* $
                    total( reform(alpha_star[i, *])*n_star[*, tm] ) & $
            endfor & $
        endfor & $
        plot_oo, tau, residual_n_star[0,*], psym=3, charsize=1.5, xstyle=1, $
            xrange=[ tau[1], tau[N_times-1] ], yrange=[1.0e-7, 1.0e-2], $
            xtitle='time !7s!3', ytitle='residual in n!u*!n and N' & $
        for i = 0, N_sizes - 1 do oplot, tau, residual_n_star[i, *], psym=3 & $
        max_residual_n_star = max(residual_n_star) & $
        print, 'max residual in equation for n_star             = ',max_residual_n_star&$

        ;Calculate residual in equation for dust abundance N.
        ;First calculate alpha_bar in units of smin_ar^2
        print, 'Checking residual in N...' & $
        alpha_bar = fltarr(N_sizes, N_sizes) & $
        for j = 0, N_sizes - 1 do begin $
            alpha_bar[0, j] = alpha_star[*, j]*((beta[j]/beta_max)^(-q_size)) $
                *((a_ring[j]/ar)^(-c_disk))/I_dust & $
        endfor & $
        ;Then calculate (dN/dt)/P_i
        dN_dt_over_P = fltarr(N_sizes, N_times) & $
        for i = 0, N_sizes - 1 do begin $
            dN_dt_over_P[i, 0] = deriv(t_sec, N[i, *])/P_i[i] & $
        endfor & $
        ;Calculate residual in the evolution equation for N
        residual_N = fltarr(N_sizes, N_times) & $
        for i = 0, N_sizes - 1 do begin $
            residual_N[i, 0] = dN_dt_over_P[i, *] - 1 & $
            for j = 0, N_sizes - 1 do begin $
                residual_N[i, 0] = residual_N[i, *] + alpha_bar[i, j]*(N[i, *]*smin_ar) $
                    *(N[j, *]*smin_ar)/(P_i[i]*Tr_sec) & $
            endfor & $
        endfor & $
        residual_N = abs(residual_N) & $
        for i = 0, N_sizes - 1 do oplot, tau, residual_N[i, *], psym=3, color=128 & $
        max_residual_N = max(residual_N) & $
        print, 'max residual in equation for N                 = ', max_residual_N & $

    endif else begin $
        residual_n_star = fltarr(N_sizes) & $
        tm = N_times - 1 & $
        for i = 0, N_sizes - 1 do $
            residual_n_star[i] = 1 - n_star_eq[i]*total(reform(alpha_star[i, *])* $
                n_star_eq)& $ 
        s = sort(beta) & $
        ymin = 1.0e-8 & $
        plot_io, beta[s], abs(residual_n_star[s])>ymin, yrange=[ymin, 1.0e-5], $
           charsize=1.5, ystyle=1, xtitle='!7b!3', $
           ytitle='fractional error in solution', psym=8 & $ 
        max_residual_n_star = max(residual_n_star) & $
        print, 'max residual in equation for n_star             = ',max_residual_n_star&$
    endelse & $
endif

;store results
save, filename='figs.dat', r, opt_depth, x, B, B0_cgs_per_ster, tau, B0_cgs_per_arc2, $
    B0_Jansky_per_arc2, s_microns, s_microns_total, N_final, N_total_final, tau, $
    A_disk_ringarea, t_col_yr, M_disk_smallestgrain, mass_gm, N_sizes, N_longitudes, $
    N_rings, beta, inputs, constants, max_residual_N, max_residual_n_star, T_a_yr, $
    T_e_yr, P_i, Tr_sec

 