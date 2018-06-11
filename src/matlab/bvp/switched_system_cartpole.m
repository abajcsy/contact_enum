function switched_system_cartpole()

    % y = (s,
    %      \theta,
    %      \dot{s},
    %      \dot{\theta},
    %      p_1,
    %      p_2,
    %      p_3,
    %      p_4,
    %      \partial s / \partial x_{n+1},
    %      \partial \theta / \partial x_{n+1},
    %      \partial \dot{s} / \partial x_{n+1},
    %      \partial \dot{\theta} / \partial x_{n+1},
    %      \partial p_1 / \partial x_{n+1},
    %      \partial p_2 / \partial x_{n+1},
    %      \partial p_3 / \partial x_{n+1},
    %      \partial p_4 / \partial x_{n+1},
    %      \lambda) TODO removing now for debugging
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
%     d = 4;
    d = 0;
    g = 9.81;
    t0 = 0;
    tfin = 8;
    
    % Initial state.
    s0 = 0;
    theta0 = 0.5;
    ds0 = 0;
    dtheta0 = 0;
    
%     t1 = 7;
    t1 = 5;
    sol = solvebvp(t1);
    
    t = evalt(sol.x, t1);
    u = evalu(sol.x, sol.y);
    
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

    function sol = solvebvp(tswitch)
        function y = guess(x)
            sfin = 1;
            thetafin = 0;
            dsfin = 0;
            dthetafin = 0;
            
%             y = zeros(18, 1);
            y = zeros(17, 1);

            if 0 <= x && x < 1
                % \tau \in [0, 1) 
                y(1) = s0 + x * (sfin - s0);
                y(2) = theta0 + x * (thetafin - theta0);
                y(3) = ds0 + x * (dsfin - ds0);
                y(4) = dtheta0 + x * (dthetafin - dtheta0);
            else
                % \tau \in [1, 2]
                y(1) = s0;
                y(2) = theta0;
                y(3) = ds0;
                y(4) = dtheta0;
            end
        end
        
        % Initialize the BVP solution. The first argument is the
        % initialization points in x, and the second argument is a function
        % that generates corresponding guesses for y
        solinit = bvpinit(linspace(0, 2, 20), @guess);

        % Solve the BVP.
        sol = bvp4c(@odefun, @bcfun, solinit);
    
        function dydx = odefun(x, y)

            s = y(1);
            theta = y(2);
            ds = y(3);
            dtheta = y(4);
            p1 = y(5);
            p2 = y(6);
            p3 = y(7);
            p4 = y(8);
            dsdswitch = y(9);
            dthetadswitch = y(10);
            ddsdswitch = y(11);
            ddthetadswitch = y(12);
            dp1dswitch = y(13);
            dp2dswitch = y(14);
            dp3dswitch = y(15);
            dp4dswitch = y(16);
            lmda1 = y(17);
%             lmda2 = y(18);
            
            if 0 <= x && x < 1
                % \tau \in [0, 1) (the not-in-contact mode)
                u1 = -(l * p3 - p4 * cos(theta))/(l * wu * (mc + mp - mp * cos(theta)^2));
                u2 = -(mc * p4 + mp * p4 - l * mp * p3 * cos(theta))/(l * mp * wu * (mc + mp - mp * cos(theta)^2));
                dydx = cartpole_odefun1(d, ddsdswitch, ddthetadswitch, dp1dswitch, dp2dswitch, dp3dswitch, dp4dswitch, ds, dtheta, dthetadswitch, g, l, mc, mp, p1, p2, p3, p4, t0, theta, tswitch, u1, u2, wp, wu);
                dydx = [dydx; 0];
%                 dydx = [dydx; 0; 0]; % \lambda does not vary over time
            else
                % \tau \in [1, 2] (the in-contact mode)
                % Note that u_1 = u_2 = 0
                dydx = cartpole_odefun2(ddsdswitch, ddthetadswitch, dp1dswitch, dp2dswitch, ds, dtheta, dthetadswitch, p1, p2, tfin, theta, tswitch, wp);
                dydx = [dydx; 0];
%                 dydx = [dydx; 0; 0]; % \lambda does not vary over time
            end
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
            dsdswitch_a = ya(9);
            dthetadswitch_a = ya(10);
            ddsdswitch_a = ya(11);
            ddthetadswitch_a = ya(12);
            dp1dswitch_a = ya(13);
            dp2dswitch_a = ya(14);
            dp3dswitch_a = ya(15);
            dp4dswitch_a = ya(16);
            lmda1_a = ya(17);
%             lmda2_a = ya(18);
            
            s_b = yb(1);
            theta_b = yb(2);
            ds_b = yb(3);
            dtheta_b = yb(4);
            p1_b = yb(5);
            p2_b = yb(6);
            p3_b = yb(7);
            p4_b = yb(8);
            dsdswitch_b = yb(9);
            dthetadswitch_b = yb(10);
            ddsdswitch_b = yb(11);
            ddthetadswitch_b = yb(12);
            dp1dswitch_b = yb(13);
            dp2dswitch_b = yb(14);
            dp3dswitch_b = yb(15);
            dp4dswitch_b = yb(16);
            lmda1_b = yb(17);
%             lmda2_b = yb(18);
            
            lmda2_a = 0;
            res = cartpole_bcfun(ddsdswitch_a,ddthetadswitch_a,dp1dswitch_b,dp2dswitch_b,dp3dswitch_b,dp4dswitch_b,ds0,ds_a,ds_b,dsdswitch_a,dtheta0,dtheta_a,dtheta_b,dthetadswitch_a,l,lmda1_a,lmda2_a,p1_b,p2_b,p3_b,p4_b,s0,s_a,s_b,sgoal,theta0,theta_a,theta_b,wg);
            res = res(1:17);
            
%             res = cartpole_bcfun(ddsdswitch_a,ddthetadswitch_a,dp1dswitch_b,dp2dswitch_b,dp3dswitch_b,dp4dswitch_b,ds0,ds_a,dsdswitch_a,dtheta0,dtheta_a,dthetadswitch_a,l,lmda_a,p1_b,p2_b,p3_b,p4_b,s0,s_a,s_b,sgoal,theta0,theta_a,theta_b,wg);
%             res = cartpole_bcfun(ddsdswitch_a,ddthetadswitch_a,dp1dswitch_b,dp2dswitch_b,dp3dswitch_b,dp4dswitch_b,ds0,ds_a,dsdswitch_a,dtheta0,dtheta_a,dthetadswitch_a,lmda_a,p1_b,p2_b,p3_b,p4_b,s0,s_a,s_b,sgoal,theta0,theta_a,wg);
%             res = [res; 0];
            %             res = cartpole_bcfun(ddsdswitch_a, ddthetadswitch_a, dp1dswitch_b, dp2dswitch_b, dp3dswitch_b, dp4dswitch_b, ds0, ds_a, dsdswitch_a, dtheta0, dtheta_a, dthetadswitch_a, l, lmda_a, lmda_b, p1_b, p2_b, p3_b, p4_b, s0, s_a, s_b, sgoal, theta0, theta_a, theta_b, wg);
%             res = res(1:16);
        end
    end

    function t = evalt(x, tswitch)
        % Converts from \tau to t using
        %   t = t_0 + (x_{n+1} - t_0) \tau, 0 \leq \tau \leq 1
        %   t = x_{n+1} + (t_f - x_{n+1})(\tau - 1), 1 \leq \tau \leq 2
        t = zeros(1, size(x, 2));
        for i = 1:size(x, 2)
            tau = x(1, i);
            
            if 0 <= tau && tau < 1
                % \tau \in [0, 1)
                t(1, i) = tswitch * tau;
            else
                % \tau \in [1, 2]
                t(1, i) = tswitch + (tfin - tswitch) * (tau - 1);
            end
        end
    end

    function u = evalu(x, y)
        % Gets u at the sampled \tau 
        u = zeros(2, size(x, 2));
        for i = 1:size(x, 2)
            tau = x(1, i);
            
            theta = y(2, i);
            p3 = y(7, i);
            p4 = y(8, i);
            
            if 0 <= tau && tau < 1
                % \tau \in [0, 1)
                u(1, i) = -(l * p3 - p4 * cos(theta))/(l * wu * (mc + mp - mp * cos(theta)^2));
                u(2, i) = -(mc * p4 + mp * p4 - l * mp * p3 * cos(theta))/(l * mp * wu * (mc + mp - mp * cos(theta)^2));
            else
                % \tau \in [1, 2]
                u(1, i) = 0;
                u(2, i) = 0;
            end
        end
    end
end