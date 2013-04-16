;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ODEINT is a Runge-Kutta integrator with adaptive stepsize control, for
;advancing y(x) from starting point x=x1 to x=x2 with fractional accuracy eps. The
;initial trial stepsize is h1, and hmin is the minumum allowed stepsize (which can
;be zero). The value of ystart is updated for x=x2. The user must supply the name of
;the procedure that calculates the derivative dydx=dy/dx from x, y, plus additional
;optional parameters that can be passed via the structure params. Adapted by Joe Hahn
;(jhahn@spacescience.org) on July 17, 2006 from Numerical Recipes in C, second edition,
;by Press et al (1992). This procedure can be called by int_sys_eqns to integrate
;a system of equations over a list of x values.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro odeint, ystart, x1, x2, eps, h1, hmin, pro_name, params=params

;INPUTS
;  ystart = scalar or vector containing the known values of y for initial value of x1
;  x1 = starting value of x
;  x2= final value of x
;  eps = fractional accuracy desired in the solution
;  h1 = initial trial stepsize
;  hmin = smallest allowed stepsize
;  pro_name = string name of proceedure that calculates dy/dx,
;       eg, pro odeint_derivs, x, y, dydx, params=params
;  params = optional structure used to pass parameters to pro_name 
;
;OUTPUTS
;  ystart = y(x) evaluated at x=x2

;constants
MAXSTP=1000000l
TINY=1.0e-30

;initialize
x=x1
h=h1
y=ystart

;take at most MAXSTP steps to perform this integration
nstp=1
finished=0
while ((nstp lt MAXSTP) and (finished eq 0)) do begin 

  ;set scaling for accuracy monitoring
  call_procedure, pro_name, x, y, dydx, params=params
  yscal = abs(y) + abs(h*dydx) + TINY

  ;decrease stepsize if next step will overshoot endpoint x2
  if ( (x+h-x2)*(x+h-x1) gt 0 ) then h=x2-x
  rkqs, y, dydx, x, h, eps, yscal, hdid, hnext, pro_name, params=params

  ;check to see if finished
  if ( (x-x2)*(x2-x1) ge 0 ) then begin
    ystart=y
    finished=1
  endif

  ;update the next step
  if (abs(hnext) lt hmin) then print,'ERROR: stepsize too small in ODEINT()'
  h=hnext
  nstp=nstp+1

endwhile
if (nstp ge MAXSTP) then print,'ERROR: too many steps in ODEINT()'

;return results
return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;RKQS is a stepsize predictor for a 5th-order Runga-Kutta integrator, with the
;stepsize chosen so that desired accuracy is achieved. The user must supply the
;procedure derivs( x, y, dydx), which calculates the derivative dydx=dy/dx. From
;Numerical Recipes in C, 2nd edition.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro rkqs, y, dydx, x, htry, eps, yscal, hdid, hnext, pro_name, params=params

;INPUTS
;  y = scalar or vector containing the known values of y at time/position x
;  dydx = derivatives dy/dx
;  x = independant coordinate
;  htry = suggested stepsize over which the solution y(x) is advanced to y(x+h)
;  eps = fractional accuracy desired in the solution
;  yscal = scale factor against which the error is scaled
;
;OUTPUTS
;  hdid = stepsize actually accomplished
;  hnext = estimated size on next step
;  x = independant time/spatial coordinate that gets updated
;  y = dependant coodinate(s), also updated by this subroutine

;some constants
safety=0.9d
pgrow=-0.2d
pshrink=-0.25d
errcon=1.89d-4

;set stepsize=initial trial value
h=htry

;loop until error threshold is satisfied
while (safety gt 0d) do begin

  ;take a step
  rkck, y, dydx, x, h, ytemp, yerr, pro_name, params=params

  ;evaluate accuracy, scaled to the required tolerance
  errmax=max(abs(yerr/yscal))/eps

  if (errmax gt 1d) then begin

    ;reduce stepsize, but by no more than a factor of 10
    htemp=safety*h*(errmax^pshrink)
    if (abs(htemp) lt abs(h*0.1d)) then htemp=h*0.1d
    h=htemp
    xnew=x+h
    if (xnew eq x) then print,'ERROR: stepsize underflow in rkqs()'

  endif else begin

   ;increase stepsize, but by no more than a factor of 5
   if (errmax gt errcon) then begin
      hnext=safety*h*(errmax^pgrow)
    endif else begin
      hnext=h*5d
    endelse
    hdid=h
    x=x+h
    y=ytemp
    break

  endelse

endwhile

;return results
return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;RKCK...Given a scalar y and its derivatives dydx evaluated for a known value of x,
;this procedure uses the 5th order Cash-Karp Runge-Kutta method to advance y over
;an interval h. The solution is returned as yout, and error estimate provided in yerr.
;The user must supply the procedure derivs( x, y, dydx), which calculates the
;derivative dydx=dy/dx. From Numerical Recipes in C, 2nd edition.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro rkck, y, dydx, x, h, yout, yerr, pro_name, params=params

;given a scalar y and its derivatives dydx evaluated for a known value of x,
;this Numerical Recipes subroutine uses the 5th order Cash-Karp 
;Runge-Kutta method to advance y over an interval h. The solution is
;returned as yout, and error estimate provided in yerr.
;The user must supply the procedure derivs( x, y, dydx), which calculates
;the derivative dydx=dy/dx
;
;INPUTS
;  y = scalar or vector containing the known values of y at time/position x
;  dydx = derivatives dy/dx
;  x = independant coordinate
;  h = small step over which the solution y(x) is advanced to yout=y(x+h)
;
;OUTPUTS
;  yout = scalar or vector containing the solution yout=y(x+h)
;  yerr = estimated error in this solution

;define some constants
A2=0.2d
A3=0.3d
A4=0.6d
A5=1d
A6=0.875d
B21=0.2d
B31=0.075d
B32=0.225d
B41=0.3d
B42=-0.9d
B43=1.2d
B51=-11d/54d
B52=2.5d
B53=-70d/27d
B54=35d/27d
B61=1631d/55296d
B62=175d/512d
B63=575d/13824d
B64=44275d/110592d
B65=253d/4096d
C1=37d/378d
C3=250d/621d
C4=125d/594d
C6=512d/1771d
DC1=C1-2825d/27648d
DC3=C3-18575d/48384d
DC4=C4-13525d/55296d
DC5=-277d/14336d
DC6=C6-0.25d

;Cash-Karp Runge-Kutta step
ytemp=y + B21*h*dydx
xtemp=x + A2*h
call_procedure, pro_name, xtemp, ytemp, ak2, params=params
ytemp=y + h*(B31*dydx + B32*ak2)
xtemp=x + A3*h
call_procedure, pro_name, xtemp, ytemp, ak3, params=params
ytemp=y + h*(B41*dydx + B42*ak2 + B43*ak3)
xtemp=x + A4*h
call_procedure, pro_name, xtemp, ytemp, ak4, params=params
ytemp=y + h*(B51*dydx + B52*ak2 + B53*ak3 + B54*ak4)
xtemp=x + A5*h
call_procedure, pro_name, xtemp, ytemp, ak5, params=params
ytemp=y + h*(B61*dydx + B62*ak2 + B63*ak3 + B64*ak4 + B65*ak5)
xtemp=x + A6*h
call_procedure, pro_name, xtemp, ytemp, ak6, params=params

;solution
yout=y + h*(C1*dydx + C3*ak3 + C4*ak4 + C6*ak6)

;error estimate = difference between 4th and 5th order methods
yerr=h*(DC1*dydx + DC3*ak3 + DC4*ak4 + DC5*ak5 + DC6*ak6)

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Perform a Runga Kutta integration of a system of equations
;y=[y0(x), y1(x), y2(x),...], and sample the output at various times x. The inputs
;are x (which is a 1D vector where yi(x) is sampled, the boundary conditions
;y0 (which is a 1D vector that are the initial yi(x[0]) values), the allowed
;fractional error, the suggested initial timestep dx_try, the minimum timestep dx_min,
;pro_name=string that provides the name of the proceedure that calculates the 
;derivative dyi/dx (a 1D vector), and the 2D output array y[j,k], where j
;indicates the quantity yi which is evaluated at x[k]. Additional parameters
;are passes on to pro_name via params, which is often a structure. Set optional
;quiet=1 to avoid printing the code's %progress.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro int_sys_eqns, x, y0, error, dx_try, dx_min, pro_name, y, params=params, quiet=quiet

;create storage for y[j, x] array
Nx = n_elements(x)
Ny = n_elements(y0)
ytype = size(y0, type=1)
y = make_array(Ny, Nx, type=ytype)
fraction = 0

;store initial conditions
y[0, 0] = y0

;loop over all x and advance y[j] from x[j-1] to x[j]
for j = 1, Nx - 1 do begin $
    ystart = y[*, j - 1] & $
    x1 = x[j - 1] & $
    x2 = x[j] & $
    odeint, ystart, x1, x2, error, dx_try, dx_min, pro_name, params=params & $
    y[0, j] = ystart & $
    if (keyword_set(quiet) eq 0) then begin $
        fraction_old = fraction & $
        fraction = round(10.0*(x2 - x[0])/(x[Nx - 1] - x[0])) & $
        if (fraction ne fraction_old) then $
            print, format = '($, A0)', strtrim(string(10*fraction),2)+'% ' & $
    endif & $
endfor
if (keyword_set(quiet) eq 0) then print, '100%'

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

