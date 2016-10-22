// napl_percolate: scilab script by walt mcnab

// a method-of-lines solution for transient 1-D gravity-driven 
// infiltration of NAPL into a partially water-saturated soil column
// water is initially present at residual saturation, but can also be
// included in an external source term

// define constants and parameters

pathway = pwd() + '\';		// local directory (to direct output)

g = 9.807; 			// gravitational acceleration, m/sec^2
u_w = 8.9e-4; 			// water viscosity, Pa*sec
u_n = 2.58e-3;			// NAPL viscosity, Pa*sec
rho_w = 1000.;			// water density, kg/m^3;
rho_n = 800.;			// NAPL density, kg/m^3;


// soil hydraulic properties
Sw_r = 0.15;
Sn_r = 0.02; 			// NAPL saturation below which its mobility drops to zero
phi = 0.35;    
k = 3.676e-13; 			// absolute permeability, m^2
n = 2.0; 			// Van Genuchten unsaturated flow shape factor
m = 1.0 - 1.0/n;


// recharge rates - placeholder
function x = q_w(t, q) 		// for water
    x = q;
endfunction
function x = q_n(t) 		// for NAPL (set to zero in this example)
    x = 0.0;
endfunction


function x = kr_w(Sw)
    // relative permeability with respect to water
    Sw = Sw.*bool2s((Sw<=1.0) & (Sw>=0.0)) + 0.0*bool2s(Sw<0.0) + 1.0*bool2s(Sw > 1.0);
    Sw_e = (Sw - Sw_r)/(1.0 - Sw_r);
    Sw_e = Sw_e. * bool2s(Sw > Sw_r);       		// prevent negative effective saturation
    x = (Sw_e.^0.5). * (1.0 - (1.0 - Sw_e.^(1.0/m)).^m).^2.0;
endfunction


function x = ambient_sat(q)
    // return saturation that allows effective hydraulic conductivity = q
    if q > k*rho_w*g/u_w then
        x = 1.0;
    else
        r = q*u_w/(k*rho_w*g);
        x0 = 0.9999;                            // initial guess
        x = fsolve(x0, list(residual_kr_w, r));
    end
endfunction
    

function x = residual_kr_w (s, r)
    // function called by fsolve to return the water saturation required to achieve kr_w = r
    x = kr_w(s) - r;
endfunction


function x = kr_n(Sn, Sw)

    // relative permeability with respect to NAPL
    Sw = Sw.*bool2s((Sw<=1.0) & (Sw>=0.0)) + 0.0*bool2s(Sw<0.0) + 1.0*bool2s(Sw > 1.0);
    Sn = Sn.*bool2s((Sn<=1.0) & (Sn>=0.0)) + 0.0*bool2s(Sn<0.0) + 1.0*bool2s(Sn > 1.0);
 
    Sw_e = (Sw - Sw_r)/(1.0 - Sw_r);
    Sw_e = Sw_e. * bool2s(Sw > Sw_r);       	// prevent negative effective saturation

    St = (Sw + Sn - Sw_r)/(1.0 - Sw_r);
    St = St. * bool2s(Sw > Sw_r);         

    x = ((St - Sw_e).^0.5). * ((1.0 - Sw_e.^(1.0/m)).^m - (1.0 - St.^(1.0/m)).^m).^2.0;
    x = x. * bool2s(Sn > Sn_r);		// reduce relative permeability to zero if S_n < residual
    
endfunction


function [ydot]=f(t,y,q, N_x,dx)

	// water volume balance
    w1 = 1/(phi*dx) * (q_w(t, q) - (k*rho_w*g/u_w)*kr_w(y(1)));
    w2 = 1/(phi*dx) * (k*rho_w*g/u_w) * (kr_w(y(1:N_x-1)) - kr_w(y(2:N_x)));
    w = [w1; w2];
    
    // NAPL volume balance
    v1 = 1/(phi*dx) * (q_n(t) - (k*rho_n*g/u_n)*kr_n(y(N_x+1), y(1)));
    v2 = 1/(phi*dx) * (k*rho_n*g/u_n) * (kr_n(y(N_x+1:$-1), y(1:N_x-1)) - kr_n(y(N_x+2:$), y(2:N_x)));
    v = [v1; v2];    

    ydot = [w; v];
    
endfunction


function xf = profile(q, L, N_x, t)	// note: user runs this function from SciLab console; must supply vaues for q, L, N_x, and t in SI units

    // ambient percolation and saturation
    S0 = ambient_sat(q);

    // column discretization 
    dx = L/N_x;
    
    // initial conditions
    for i = 1:N_x,
        x(i) = (i-0.5)*dx;
        y0(i) = S0;      		// water saturation
        y0(N_x + i) = 0.0       	// NAPL saturation = zero everywhere
    end
    
    // initial conditions for NAPL spill
    for i = 1:4,
        y0(N_x + i) = 0.58       	// NAPL saturation = zero everywhere
    end    

    t0 = 0.0 							// simulation starting time
    %ODEOPTIONS=[1,0,0,%inf,0,2,4000,12,5,0,-1,-1];             // set maxsteps for ODE solver = 4000
    y = ode(y0,t0,t,list(f, q, N_x, dx));
    
    // plot solution
    xtitle('Saturation', 'z', 'S');
    plot2d(x, [y(1:N_x) y(N_x+1:$)], leg="Water@NAPL");
    
    // write to file
    a_out = [x y(1:N_x) y(N_x+1:$)];
    file_id = pathway + 'saturation.txt'
    deletefile(file_id);  
    write(file_id, a_out, '(e11.5, 2x, e11.5, 2x, e11.5)');
   
   // return location of NAPL front
   S_min = 0.001;                   // threshold NAPL saturation for quantifying location of front
   loc = find(y(N_x+1:$) > S_min);
   xf = (max(loc)-0.5)*dx;
    
endfunction

