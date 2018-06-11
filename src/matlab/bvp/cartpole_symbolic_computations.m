
% fprintf('## Computing dynamics...\n');
% compute_cartpole_dynamics()
% fprintf('## ...done.\n');

% ----- Begin generate ----- %
% fprintf('## Computing subsystem 1 ODE...\n');
% compute_ode1()
% fprintf('## ...done.\n');

% fprintf('## Computing subsystem 2 ODE...\n');
% compute_ode2()
% fprintf('## ...done.\n');

fprintf('## Computing boundary conditions...\n');
compute_bc()
fprintf('## ...done.\n');

% fprintf('## Computing partial J_1 wrt switching time...\n');
% compute_dJ();
% fprintf('## ...done.\n');
% ----- End generate ----- %

% ----- Begin tests ----- %
% fprintf('## Running test on subsystem 1 ODE...\n');
% test_ode1()
% fprintf('## ...done.\n');
% ----- End tests ----- %

function compute_cartpole_dynamics()
    syms mc mp l g d u1 u2 s theta ds dtheta dds ddtheta;
    assume([mc mp l g d u1 u2 s theta ds dtheta dds ddtheta], 'real');

    % Solve for the dynamics in a simpler form. 
    % ddq = [(1/(mc + mp)) * ((d + u(1)) - mp * l * (ddq(2) * cos(q(2)) - dq(2)^2 * sin(q(2))));
    %        u(2)/(l * mp) + (g/l) * sin(q(2)) - (ddq(1)/l) * cos(q(2))];

    eqns = [dds == (1/(mc + mp)) * ((d + u1) - mp * l * (ddtheta * cos(theta) - dtheta^2 * sin(theta))), ddtheta == u2/(l * mp) + (g/l) * sin(theta) - (dds/l) * cos(theta)];
    vars = [dds, ddtheta];
    soln = solve(eqns, vars);
    soln.dds
    soln.ddtheta
end

function compute_ode1()
    % Computes the ODE functions symbolically. 
    % Not-in-contact subsystem. 

    syms mc mp l g d u1 u2 s theta ds dtheta dds ddtheta tswitch t0 p1 p2 p3 p4 wu wp;
    assume([mc mp l g d u1 u2 s theta ds dtheta dds ddtheta tswitch t0 p1 p2 p3 p4 wu wp], 'real');
    
    % x = (s, \theta, \dot{s}, \dot{\theta})^\top
    % u = (u_1, u_2)^\top
    
    % Cartpole dynamics computed above, where \dot{x} = f(x, u)
    f = [ds;
         dtheta;
         (l*mp*sin(theta)*dtheta^2 + d + u1 - u2*cos(theta) - g*mp*cos(theta)*sin(theta))/(mc + mp - mp*cos(theta)^2);
         (mc*u2 + mp*u2 - d*mp*cos(theta) - mp*u1*cos(theta) + g*mp^2*sin(theta) + g*mc*mp*sin(theta) - dtheta^2*l*mp^2*cos(theta)*sin(theta))/(l*mp*(mc + mp - mp*cos(theta)^2))];
    ftilde = (tswitch - t0) * f;
    
    % The running cost. 
    L = 0.5 * (wu * (u1^2 + u2^2) + wp * theta^2);
    Ltilde = (tswitch - t0) * L;
    
    % Compute \partial x \partial \tau
    dxdtau = ftilde;
    
    % Compute \partial p / \partial \tau with p = (p_1, p_2, p_3, p_4)^\top
    p = [p1; p2; p3; p4];
    
    dftildedx = [diff(ftilde, s), diff(ftilde, theta), diff(ftilde, ds), diff(ftilde, dtheta)];
             
    dLtildedx = [diff(Ltilde, s), diff(Ltilde, theta), diff(Ltilde, ds), diff(Ltilde, dtheta)];
             
    dpdtau = -dftildedx' * p - dLtildedx';
             
    % Solve for u_1 and u_2 in terms of p_1, p_2, p_3, and p_4
    dftildedu = [diff(ftilde, u1), diff(ftilde, u2)];
             
    dLtildedu = [diff(Ltilde, u1), diff(Ltilde, u2)];
             
    res = dftildedu' * p + dLtildedu';
    eqns = [0 == res(1), 0 == res(2)];
    vars = [u1, u2];
    soln = solve(eqns, vars);
    fprintf('u1:\n');
    disp(soln.u1);
    fprintf('u2:\n');
    disp(soln.u2);

    syms dsdswitch dthetadswitch ddsdswitch ddthetadswitch dp1dswitch dp2dswitch dp3dswitch dp4dswitch
    assume([dsdswitch dthetadswitch ddsdswitch ddthetadswitch dp1dswitch dp2dswitch dp3dswitch dp4dswitch], 'real');
    
    % Taken from the above computation of u_1 and u_2
    du1dswitch = -((dp3dswitch*(t0 - tswitch))/(mc + mp - mp*cos(theta)^2) - (dp4dswitch*cos(theta)*(t0 - tswitch))/(l*(mc + mp - mp*cos(theta)^2)))/(wu*(t0 - tswitch));
    du2dswitch = ((dp3dswitch*cos(theta)*(t0 - tswitch))/(mc + mp - mp*cos(theta)^2) - (dp4dswitch*(mc + mp)*(t0 - tswitch))/(l*mp*(mc + mp - mp*cos(theta)^2)))/(wu*(t0 - tswitch));
    
    dxdswitch = [dsdswitch;
                 dthetadswitch;
                 ddsdswitch;
                 ddthetadswitch];
             
    dudswitch = [du1dswitch;
                 du2dswitch];
             
    dpdswitch = [dp1dswitch;
                 dp2dswitch;
                 dp3dswitch;
                 dp4dswitch];
             
    dfdx = [diff(f, s), diff(f, theta), diff(f, ds), diff(f, dtheta)];
    dfdu = [diff(f, u1), diff(f, u2)];

    % Compute \partial / \partial \tau ( \partial x / \partial x_{n+1})
    tau_dxdswitch = f + (tswitch - t0) * (dfdx * dxdswitch + dfdu * dudswitch);
    
    dLdx = [diff(L, s), diff(L, theta), diff(L, ds), diff(L, dtheta)];
    
    dLdu = [diff(L, u1), diff(L, u2)];
    
    term1 = dfdx' * dpdswitch;
    
    % For convenience define x, u vectors
    x = [s; theta; ds; dtheta];
    u = [u1; u2];
    n = 4;
    m = 2;
    
    term2 = sym(zeros(4, 1));
    
    for j2 = 1:n
        for j1 = 1:n
            for j3 = 1:n
                term2(j2) = term2(j2) + p(j1) * diff(diff(f(j1), x(j3)), x(j2)) * dxdswitch(j3);
            end
        end
    end
    
    term3 = sym(zeros(4, 1));
    
    for j2 = 1:n
        for j1 = 1:n
            for j3 = 1:m
                term3(j2) = term3(j2) + p(j1) * diff(diff(f(j1), u(j3)), x(j2)) * dudswitch(j3);
            end
        end
    end
        
    term4 = [diff(dLdx, s); 
             diff(dLdx, theta); 
             diff(dLdx, ds); 
             diff(dLdx, dtheta)] * dxdswitch;
 
    % This term should be 0.
    term5 = [diff(dLdu, s);
             diff(dLdu, theta);
             diff(dLdu, ds);
             diff(dLdu, dtheta)] * dudswitch;
    
    terms = term1 + term2 + term3 + term4 + term5;
    
    tau_dpdswitch = -dfdx' * p - dLdx' - (tswitch - t0) * terms;
    
    % Output the results.
    dydx = [dxdtau;
            dpdtau;
            tau_dxdswitch;
            tau_dpdswitch];
    matlabFunction(dydx, 'File', 'cartpole_odefun1');
end

function compute_ode2()
    % Computes the ODE functions symbolically. 
    % In-contact subsystem. 

    syms mc mp l g d u1 u2 s theta ds dtheta dds ddtheta tswitch tfin p1 p2 p3 p4 wu wp;
    assume([mc mp l g d u1 u2 s theta ds dtheta dds ddtheta tswitch tfin p1 p2 p3 p4 wu wp], 'real');
    
    % x = (s, \theta, \dot{s}, \dot{\theta})^\top
    % u = (u_1, u_2)^\top
    
    % Cartpole dynamics when in contact with the wall
    f = [ds;
         dtheta;
         0;
         0];
    ftilde = (tfin - tswitch) * f;
    
    % The running cost. 
    L = 0.5 * (wu * (u1^2 + u2^2) + wp * theta^2);
    Ltilde = (tfin - tswitch) * L;
    
    % Compute \partial x \partial \tau
    dxdtau = ftilde;
    
    % Compute \partial p / \partial \tau with p = (p_1, p_2, p_3, p_4)^\top
    p = [p1; p2; p3; p4];
    
    dftildedx = [diff(ftilde, s), diff(ftilde, theta), diff(ftilde, ds), diff(ftilde, dtheta)];
             
    dLtildedx = [diff(Ltilde, s), diff(Ltilde, theta), diff(Ltilde, ds), diff(Ltilde, dtheta)];
             
    dpdtau = -dftildedx' * p - dLtildedx';
             
    % Solve for u_1 and u_2 in terms of p_1, p_2, p_3, and p_4
    dftildedu = [diff(ftilde, u1), diff(ftilde, u2)];
             
    dLtildedu = [diff(Ltilde, u1), diff(Ltilde, u2)];
             
    res = dftildedu' * p + dLtildedu';
    eqns = [0 == res(1), 0 == res(2)];
    vars = [u1, u2];
    soln = solve(eqns, vars);
    fprintf('u1:\n');
    disp(soln.u1);
    fprintf('u2:\n');
    disp(soln.u2);

    syms dsdswitch dthetadswitch ddsdswitch ddthetadswitch dp1dswitch dp2dswitch dp3dswitch dp4dswitch
    assume([dsdswitch dthetadswitch ddsdswitch ddthetadswitch dp1dswitch dp2dswitch dp3dswitch dp4dswitch], 'real');
    
    % Taken from the above computation of u_1 and u_2
    du1dswitch = 0;
    du2dswitch = 0;
    
    dxdswitch = [dsdswitch;
                 dthetadswitch;
                 ddsdswitch;
                 ddthetadswitch];
             
    dudswitch = [du1dswitch;
                 du2dswitch];
             
    dpdswitch = [dp1dswitch;
                 dp2dswitch;
                 dp3dswitch;
                 dp4dswitch];
             
    dfdx = [diff(f, s), diff(f, theta), diff(f, ds), diff(f, dtheta)];
    dfdu = [diff(f, u1), diff(f, u2)];

    % Compute \partial / \partial \tau ( \partial x / \partial x_{n+1})
    tau_dxdswitch = -f + (tfin - tswitch) * (dfdx * dxdswitch + dfdu * dudswitch);
    
    dLdx = [diff(L, s), diff(L, theta), diff(L, ds), diff(L, dtheta)];
    
    dLdu = [diff(L, u1), diff(L, u2)];
    
    term1 = dfdx' * dpdswitch;
    
    % For convenience define x, u vectors
    x = [s; theta; ds; dtheta];
    u = [u1; u2];
    n = 4;
    m = 2;
    
    term2 = sym(zeros(4, 1));
    
    for j2 = 1:n
        for j1 = 1:n
            for j3 = 1:n
                term2(j2) = term2(j2) + p(j1) * diff(diff(f(j1), x(j3)), x(j2)) * dxdswitch(j3);
            end
        end
    end
    
    term3 = sym(zeros(4, 1));
    
    for j2 = 1:n
        for j1 = 1:n
            for j3 = 1:m
                term3(j2) = term3(j2) + p(j1) * diff(diff(f(j1), u(j3)), x(j2)) * dudswitch(j3);
            end
        end
    end
        
    term4 = [diff(dLdx, s); 
             diff(dLdx, theta); 
             diff(dLdx, ds); 
             diff(dLdx, dtheta)] * dxdswitch;
 
    term5 = [diff(dLdu, s);
             diff(dLdu, theta);
             diff(dLdu, ds);
             diff(dLdu, dtheta)] * dudswitch;
    
    terms = term1 + term2 + term3 + term4 + term5;
    
    tau_dpdswitch = dfdx' * p + dLdx' - (tfin - tswitch) * terms;
    
    % Output the results.
    dydx = [dxdtau;
            dpdtau;
            tau_dxdswitch;
            tau_dpdswitch]
    matlabFunction(dydx, 'File', 'cartpole_odefun2');
end

function compute_bc()
    % Computes the boundary conditions symbolically.
    
    % Initial state x_0 = (s_0, \theta_0, \dot{s}_0, \dot{\theta}_0)^\top 
    % Note that x(0) = x(a) = x_0
    
    syms s_a theta_a ds_a dtheta_a s0 theta0 ds0 dtheta0 s_b theta_b ds_b dtheta_b dsdswitch_a dthetadswitch_a ddsdswitch_a ddthetadswitch_a p1_b p2_b p3_b p4_b dp1dswitch_b dp2dswitch_b dp3dswitch_b dp4dswitch_b wg sgoal swall l lmda1_a lmda2_a lmda1_b lmda2_b;
    assume([s_a theta_a ds_a dtheta_a s0 theta0 ds0 dtheta0 s_b theta_b ds_b dtheta_b dsdswitch_a dthetadswitch_a ddsdswitch_a ddthetadswitch_a p1_b p2_b p3_b p4_b dp1dswitch_b dp2dswitch_b dp3dswitch_b dp4dswitch_b wg sgoal swall l lmda1_a lmda2_a lmda1_b lmda2_b], 'real');
    
    % The final cost.
    psi = wg * (s_b - sgoal)^2;
    dpsidx = [diff(psi, s_b), diff(psi, theta_b), diff(psi, ds_b), diff(psi, dtheta_b)];
    ddpsidx = [diff(dpsidx, s_b);
               diff(dpsidx, theta_b);
               diff(dpsidx, ds_b);
               diff(dpsidx, dtheta_b)];
    
    % The final state constraint.
%     phi = swall - (l * sin(theta_b) + s_b);
    phi = [swall - (l * sin(theta_b) + s_b);
           ds_b^2 + dtheta_b^2];       
    dphidx = [diff(phi, s_b), diff(phi, theta_b), diff(phi, ds_b), diff(phi, dtheta_b)];
    ddphidx = [diff(dphidx, s_b);
               diff(dphidx, theta_b);
               diff(dphidx, ds_b);
               diff(dphidx, dtheta_b)];
           
    % Compute \lambda
    p = [p1_b; p2_b; p3_b; p4_b];
    lmda = [lmda1_a; lmda2_a];
%     eqns = p == dpsidx' + dphidx' * lmda
%     soln = solve(eqns, lmda, 'ReturnConditions', true)
    % TODO Why does the above fail to find a solution?
    lmda_bc = [lmda1_a - ((2 * wg * (s_b - sgoal) - (p1_b + p2_b))/(1 + l * cos(theta_b)));
               0];
%                lmda2_a - ((p3_b + p4_b)/(2 * (ds_b + dtheta_b)))];
    
    p_bc = p - (dpsidx' + dphidx' * lmda);
    
    dpdswitch = [dp1dswitch_b; dp2dswitch_b; dp3dswitch_b; dp4dswitch_b];
    dpdswitch_bc = dpdswitch - (dpsidx' + dphidx' * lmda);
            
    bc = [s_a - s0;
          theta_a - theta0;
          ds_a - ds0
          dtheta_a - dtheta0;
          p_bc;
          dsdswitch_a;
          dthetadswitch_a;
          ddsdswitch_a;
          ddthetadswitch_a;
          dpdswitch_bc;
          lmda_bc]; 
    matlabFunction(bc, 'File', 'cartpole_bcfun');
end

function compute_dJ()
    syms wg sgoal s_b theta_b ds_b dtheta_b dsdswitch_b dthetadswitch_b ddsdswitch_b ddthetadswitch_b;
    assume([wg sgoal s_b theta_b ds_b dtheta_b dsdswitch_b dthetadswitch_b ddsdswitch_b ddthetadswitch_b], 'real');
    
    % Compute first term. 
    psi = wg * (s_b - sgoal)^2;
    dpsidx = [diff(psi, s_b), diff(psi, theta_b), diff(psi, ds_b), diff(psi, dtheta_b)];
    dxdswitch_b = [dsdswithc_b; dthetadswitch_b; ddsdswitch_b; ddthetadswitch_b];
    disp(dpsidx * dxdswitch_b);
    
    % Taken from the above computation of u_1 and u_2
    dudswitch_1 = [-((dp3dswitch*(t0 - tswitch))/(mc + mp - mp*cos(theta)^2) - (dp4dswitch*cos(theta)*(t0 - tswitch))/(l*(mc + mp - mp*cos(theta)^2)))/(wu*(t0 - tswitch));
                   ((dp3dswitch*cos(theta)*(t0 - tswitch))/(mc + mp - mp*cos(theta)^2) - (dp4dswitch*(mc + mp)*(t0 - tswitch))/(l*mp*(mc + mp - mp*cos(theta)^2)))/(wu*(t0 - tswitch))];
    dudswitch_2 = [0; 0];
    
    % Compute the term inside the first integrand.
    % TODO finish this
end

function test_ode1()    
    % y = (x_1, x_2, x_3, x_4, p_1, p_2, p_3, p_4)
    % x = \tau
    
    % Cost function weights.
    wp = 5;
    wu = 1;
    wg = 30;
    
    % Goal position.
    sgoal = 1;
    
    % Parameters. 
    mp = 1;
    mc = 10;
    l = 0.5;
    d = 0;
    g = 9.81;
    t0 = 0;
    tswitch = 8;
    
    % Initial state.
    s0 = 0;
%     s0 = 0.25;
    theta0 = 0;
%     theta0 = pi/8;
%     theta0 = pi/4;
    ds0 = 0;
    dtheta0 = 0;
    
    function y = guess(x)
        sfin = 1;
        thetafin = 0;
        dsfin = 0;
        dthetafin = 0;
        
        y = [s0 + x * (sfin - s0);
             theta0 + x * (thetafin - theta0);
             ds0 + x * (dsfin - ds0);
             dtheta0 + x * (dthetafin - dtheta0);
             0;
             0;
             0;
             0];
    end
    
%     solinit = bvpinit([0, 1], @guess);
    solinit = bvpinit(linspace(0, 1, 10), @guess);
    sol = bvp4c(@odefun, @bcfun, solinit);
    
    t = evalt(sol.x);
    u = evalu(sol.y);
    
    % Plot the state and controls over time.
    figure(1);
    
    subplot(2, 2, 1);
    plot(t, sol.y(1, :), 'b-x');
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$s$', 'Interpreter', 'latex');

    subplot(2, 2, 2);
    plot(t, sol.y(2, :), 'b-x');
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$\theta$', 'Interpreter', 'latex');
    
    subplot(2, 2, 3);
    plot(t, u(1, :), 'b-x');
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$u_1$', 'Interpreter', 'latex');
    
    subplot(2, 2, 4);
    plot(t, u(2, :), 'b-x');
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$u_2$', 'Interpreter', 'latex');
    
    % Visualize the cartpole trajectory over time. 
    figure(2);
    q = [s0; theta0; ds0; dtheta0];
    cart_width = 0.6;
    cart_height = 0.3;
    cart_center = [0; 0.25];
    
    axis([-2.5, 2.5, -1, 2.5]);
    
    lc1 = line([q(1) + cart_center(1) - 0.5 * cart_width, q(1) + cart_center(1) - 0.5 * cart_width], [cart_center(2) - 0.5 * cart_height, cart_center(2) + 0.5 * cart_height]);
    lc2 = line([q(1) + cart_center(1) - 0.5 * cart_width, q(1) + cart_center(1) + 0.5 * cart_width], [cart_center(2) + 0.5 * cart_height, cart_center(2) + 0.5 * cart_height]);
    lc3 = line([q(1) + cart_center(1) + 0.5 * cart_width, q(1) + cart_center(1) + 0.5 * cart_width], [cart_center(2) + 0.5 * cart_height, cart_center(2) - 0.5 * cart_height]);
    lc4 = line([q(1) + cart_center(1) + 0.5 * cart_width, q(1) + cart_center(1) - 0.5 * cart_width], [cart_center(2) - 0.5 * cart_height, cart_center(2) - 0.5 * cart_height]);
    
    lp1 = line([q(1) + cart_center(1), q(1) + cart_center(1) + l * sin(q(2))], [cart_center(2) + 0.5 * cart_height, cart_center(2) + 0.5 * cart_height + l * cos(q(2))]);

    % TODO Compute interpolated trajectory with fixed timestep (e.g. 0.01)
    % for smoother visualization. 
    
    last_t = 0;
    for k = 1:size(sol.y, 2)
        q = sol.y(1:4, k);
        dt = t(k) - last_t;
        last_t = t(k);
        
        % Update the visualization of the cartpole.
        set(lc1, 'Xdata', [q(1) + cart_center(1) - 0.5 * cart_width, q(1) + cart_center(1) - 0.5 * cart_width]);
        set(lc2, 'Xdata', [q(1) + cart_center(1) - 0.5 * cart_width, q(1) + cart_center(1) + 0.5 * cart_width]);
        set(lc3, 'Xdata', [q(1) + cart_center(1) + 0.5 * cart_width, q(1) + cart_center(1) + 0.5 * cart_width]);
        set(lc4, 'Xdata', [q(1) + cart_center(1) + 0.5 * cart_width, q(1) + cart_center(1) - 0.5 * cart_width]);
        set(lp1, 'Xdata', [q(1) + cart_center(1), q(1) + cart_center(1) + l * sin(q(2))]);
        set(lp1, 'Ydata', [cart_center(2) + 0.5 * cart_height, cart_center(2) + 0.5 * cart_height + l * cos(q(2))]);
        pause(dt);
    end    
    
    function dydx = odefun(x, y) 
        s = y(1);
        theta = y(2);
        ds = y(3);
        dtheta = y(4);
        p1 = y(5);
        p2 = y(6);
        p3 = y(7);
        p4 = y(8);
        
        u1 = -(l * p3 - p4 * cos(theta))/(l * wu * (mc + mp - mp * cos(theta)^2));
        u2 = -(mc * p4 + mp * p4 - l * mp * p3 * cos(theta))/(l * mp * wu * (mc + mp - mp * cos(theta)^2));        
        dydx = cartpole_odefun1(d, 0, 0, 0, 0, 0, 0, ds, dtheta, 0, g, l, mc, mp, p1, p2, p3, p4, t0, theta, tswitch, u1, u2, wp, wu);
        dydx = dydx(1:8);
    end

    function res = bcfun(ya, yb)
        s_a = ya(1);
        theta_a = ya(2);
        ds_a = ya(3);
        dtheta_a = ya(4);
        p1_a = ya(5);
        p2_a = ya(6);
        p3_a = ya(7);
        p4_a = ya(8);

        s_b = yb(1);
        theta_b = yb(2);
        ds_b = yb(3);
        dtheta_b = yb(4);
        p1_b = yb(5);
        p2_b = yb(6);
        p3_b = yb(7);
        p4_b = yb(8);

        res = cartpole_bcfun(0, 0, 0, 0, 0, 0, ds0, ds_a, 0, dtheta0, dtheta_a, 0, l, 0, 0, p1_b, p2_b, p3_b, p4_b, s0, s_a, s_b, sgoal, theta0, theta_a, theta_b, wg);
        res = res(1:8);
    end

    function t = evalt(x)
        % Converts from \tau to t using
        %   t = t_0 + (x_{n+1} - t_0) \tau
        t = zeros(1, size(x, 2));
        for i = 1:size(x, 2)
            tau = x(1, i);
            t(1, i) = tswitch * tau;
        end
    end

    function u = evalu(y)
        % Gets u at the sampled points
        u = zeros(2, size(y, 2));
        for i = 1:size(y, 2)

            theta = y(2, i);
            p3 = y(7, i);
            p4 = y(8, i);
            
            u(1, i) = -(l * p3 - p4 * cos(theta))/(l * wu * (mc + mp - mp * cos(theta)^2));
            u(2, i) = -(mc * p4 + mp * p4 - l * mp * p3 * cos(theta))/(l * mp * wu * (mc + mp - mp * cos(theta)^2));
        end
    end
end