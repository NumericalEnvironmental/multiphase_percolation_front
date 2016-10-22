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


function soil = soil_sample()
    // return a random soil column property configuration (represented as a structure)
    k = grand(1,1,'unf', 1e-13, 1e-12);
    phi = grand(1,1,'unf', 0.2, 0.35);
    n = grand(1,1,'unf', 1.5, 2.5);
    m = 1.0 - n.^(-1);
    Sw_r = grand(1,1,'unf', 0.075, 0.2);
    Sn_r = grand(1,1,'unf', 0.0, 0.035);    
    soil = struct('k', k, 'phi', phi, 'n', n, 'm', m, 'Sw_r', Sw_r, 'Sn_r', Sn_r);
endfunction


// recharge rates - placeholder
function x = q_w(t, q) 		// for water
    x = q;
endfunction
function x = q_n(t) 		// for NAPL; not used in this example
    x = 0.0;
endfunction


function x = kr_w(Sw, soil)
    // relative permeability with respect to water
    Sw = Sw.*bool2s((Sw<=1.0) & (Sw>=0.0)) + 0.0*bool2s(Sw<0.0) + 1.0*bool2s(Sw > 1.0);
    Sw_e = (Sw - soil.Sw_r)/(1.0 - soil.Sw_r);
    Sw_e = Sw_e. * bool2s(Sw > soil.Sw_r);       		// prevent negative effective saturation
    x = (Sw_e.^0.5). * (1.0 - (1.0 - Sw_e.^(1.0/soil.m)).^soil.m).^2.0;
endfunction


function x = ambient_sat(q, soil)
    // return saturation that allows effective hydraulic conductivity = q
    if q > soil.k*rho_w*g/u_w then
        x = 1.0;
    else
        r = q*u_w/(soil.k*rho_w*g);
        x0 = 0.9999;                            // initial guess
        x = fsolve(x0, list(residual_kr_w, r, soil));
    end
endfunction
    

function x = residual_kr_w (s, r, soil)
    // function called by fsolve to return the water saturation required to achieve kr_w = r
    x = kr_w(s, soil) - r;
endfunction


function x = kr_n(Sn, Sw, soil)

    // relative permeability with respect to NAPL
    Sw = Sw.*bool2s((Sw<=1.0) & (Sw>=0.0)) + 0.0*bool2s(Sw<0.0) + 1.0*bool2s(Sw > 1.0);
    Sn = Sn.*bool2s((Sn<=1.0) & (Sn>=0.0)) + 0.0*bool2s(Sn<0.0) + 1.0*bool2s(Sn > 1.0);
 
    Sw_e = (Sw - soil.Sw_r)/(1.0 - soil.Sw_r);
    Sw_e = Sw_e. * bool2s(Sw > soil.Sw_r);       	// prevent negative effective saturation

    St = (Sw + Sn - soil.Sw_r)/(1.0 - soil.Sw_r);
    St = St. * bool2s(Sw > soil.Sw_r);         

    x = ((St - Sw_e).^0.5). * ((1.0 - Sw_e.^(1.0/soil.m)).^soil.m - (1.0 - St.^(1.0/soil.m)).^soil.m).^2.0;
    x = x. * bool2s(Sn > soil.Sn_r);		// reduce relative permeability to zero if S_n < residual
    
endfunction


function [ydot]=f(t, y, q, soil, N_x, dx)

    // water volume balance
    w1 = 1/(soil.phi*dx) * (q_w(t, q) - (soil.k*rho_w*g/u_w)*kr_w(y(1)));
    w2 = 1/(soil.phi*dx) * (soil.k*rho_w*g/u_w) * (kr_w(y(1:N_x-1)) - kr_w(y(2:N_x)));
    w = [w1; w2];
    
    // NAPL volume balance
    v1 = 1/(soil.phi*dx) * (q_n(t) - (soil.k*rho_n*g/u_n)*kr_n(y(N_x+1), y(1)));
    v2 = 1/(soil.phi*dx) * (soil.k*rho_n*g/u_n) * (kr_n(y(N_x+1:$-1), y(1:N_x-1), soil) - kr_n(y(N_x+2:$), y(2:N_x), soil));
    v = [v1; v2];    

    ydot = [w; v];
    
endfunction


function monte_carlo(num_trials) 						// user call this function from SciLab console, supplying number of trials
    
    // problem definition - column length (e.g., 100 ft), number of cells, end time (e.g., 50 years)
    L = 100.0/3.281;
    N_x = 100;
    t = 50*365.25*86400;
    
    q = (grand(num_trials,1,'unf',0,5) * 2.54/100)/(86400*365.25);             // q = percolation, converted from in/yr (i.e., 0-5)
    for i = 1:num_trials
        soil(i) = soil_sample();
    end
    
    for i = 1:num_trials
        disp(i);        
        d(i) = profile(q(i), soil(i), L, N_x, t);
        disp([soil(i).k d(i)]);
    end
  
    // write to output file
    for i = 1:num_trials
        k(i) = soil(i).k;
        phi(i) = soil(i).phi;        
        n(i) = soil(i).n;        
        Sw_r(i) = soil(i).Sw_r;         
        Sn_r(i) = soil(i).Sn_r;         
    end
    a_out = [q k phi n Sw_r Sn_r d];
    file_id = pathway + 'mc_log.txt'
    deletefile(file_id);  
    write(file_id, a_out, '(e11.5, 2x, e11.5, 2x, e11.5, 2x, e11.5, 2x, e11.5, 2x, e11.5, 2x, e11.5)');  
    
endfunction


function xf = profile(q, soil, L, N_x, t)

    // ambient percolation and saturation
    Sw_0 = ambient_sat(q, soil);

    // column discretization 
    dx = L/N_x;
    
    // initial column-long conditions
    for i = 1:N_x,
        x(i) = (i-0.5)*dx;
        y0(i) = Sw_0;      		// water saturation
        y0(N_x + i) = 0.0       // NAPL saturation = zero everywhere
    end

    // initial NAPL saturation
    Sn_0 = 1.0 - (Sw_0 + 0.01);
    spill_cells = int(2.0/3.0 * 0.04 * N_x);
    for i = 1:spill_cells,
        y0(N_x + i) = Sn_0;
    end    

    t0 = 0.0 					// simulation starting time
    %ODEOPTIONS=[1,0,0,%inf,0,2,10000,12,5,0,-1,-1];             // set maxsteps for ODE solver = 1000
    y = ode(y0,t0,t,list(f, q, soil, N_x, dx));
    
   // return location of NAPL front
   S_min = 0.001;                   // threshold NAPL saturation for quantifying location of front
   loc = find(y(N_x+1:$) > S_min);
   xf = (max(loc)-0.5)*dx;
    
endfunction

