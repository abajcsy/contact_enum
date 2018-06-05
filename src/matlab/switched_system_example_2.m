function switched_system_example_2()

    % y = (x_1, 
    %      x_2, 
    %      p_1, 
    %      p_2, 
    %      \partial x_1 / \partial x_{n+1},
    %      \partial x_2 / \partial x_{n+1},
    %      \partial p_1 / \partial x_{n+1},
    %      \partial p_2 / \partial x_{n+1})
    % x = \tau
    
    % The switching time.
%     t1 = 0.1897;
    t1 = 1;
    
    alpha = 0.001;
    
    % Use a gradient projection method to find optimal switching time.
    for i = 1:100
        % Solve the BVP ODE.
        sol = solvebvp(t1);

        dJ = evaldJ(sol.x, sol.y, t1);
        J = evalJ(sol.x, sol.y, t1);

        fprintf('t1 = %f, J = %f, dJ = %f\n', t1, J, dJ);
        
        % TODO Add in gradient projection.
        % TODO Add in line search to choose alpha.
        
        t1 = t1 - alpha * dJ;
    end

    % Plot the results. 
    figure(1);
    plot(sol.y(1, :), sol.y(2, :), 'b-x');
    xlabel('$x_1$', 'Interpreter', 'latex');
    ylabel('$x_2$', 'Interpreter', 'latex');
    
    figure(2);
    u = evalu(sol.x, sol.y);
    t = evalt(sol.x, t1);
    plot(t(1, :), u(1, :));
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$u$', 'Interpreter', 'latex');
    
    function sol = solvebvp(tswitch)
        % Initialize the BVP solution. The first argument is the bound on \tau,
        % and the second argument is the initial guess
        solinit = bvpinit([0, 2], zeros(8, 1));

        % Solve the BVP.
        sol = bvp4c(@odefun, @bcfun, solinit);
    
        function dydx = odefun(x, y)

            x1 = y(1);
            x2 = y(2);
            p1 = y(3);
            p2 = y(4);
            dx1dswitch = y(5);
            dx2dswitch = y(6);
            dp1dswitch = y(7);
            dp2dswitch = y(8);

            if 0 <= x && x < 1
                % \tau \in [0, 1)
                u = -p1 - p2;
                dydx = [tswitch * (0.6 * x1 + 1.2 * x2 + u);            % x1
                        tswitch * (-0.8 * x1 + 3.4 * x2 + u);           % x2
                        -tswitch * (0.6 * p1 - 0.8 * p2);               % p1
                        -tswitch * (1.2 * p1 + 3.4 * p2 + x2 - 2);      % p2
                        0.6 * x1 + 1.2 * x2 + u + tswitch * (0.6 * dx1dswitch + 1.2 * dx2dswitch - dp1dswitch - dp2dswitch);    % dx1dswitch
                        -0.8 * x1 + 3.4 * x2 + u + tswitch * (-0.8 * dx1dswitch + 3.4 * dx2dswitch - dp1dswitch - dp2dswitch);  % dx2dswitch
                        -(0.6 * p1 - 0.8 * p2) - tswitch * (0.6 * dp1dswitch - 0.8 * dp2dswitch);                               % dp1dswitch
                        -(1.2 * p1 + 3.4 * p2) - (x2 - 2) - tswitch * (1.2 * dp1dswitch + 3.4 * dp2dswitch + dx2dswitch)];      % dp2dswitch
            else
                % \tau \in [1, 2]
                u = p2 - 2 * p1;
                dydx = [(2 - tswitch) * (4 * x1 + 3 * x2 + 2 * u);      % x1
                        (2 - tswitch) * (-x1 - u);                      % x2
                        -(2 - tswitch) * (4 * p1 - p2);                 % p1
                        -(2 - tswitch) * (3 * p1 + (x2 - 2));           % p2
                        -(4 * x1 + 3 * x2 + 2 * u) + (2 - tswitch) * (4 * dx1dswitch + 3 * dx2dswitch + 2 * (dp2dswitch - 2 * dp1dswitch));  % dx1dswitch
                        -(-x1 - u) + (2 - tswitch) * (-dx1dswitch - (dp2dswitch - 2 * dp1dswitch));                                          % dx2dswitch
                        4 * p1 - p2 - (2 - tswitch) * (4 * dp1dswitch - dp2dswitch);                                                         % dp1dswitch   
                        3 * p1 + (x2 - 2) - (2 - tswitch) * (3 * dp1dswitch + dx2dswitch)];                                                  % dp2dswitch
            end
        end

        function res = bcfun(ya, yb)

            x1_a = ya(1);
            x2_a = ya(2);
            dx1dswitch_a = ya(5);
            dx2dswitch_a = ya(6);

            x1_b = yb(1);
            x2_b = yb(2);
            p1_b = yb(3);
            p2_b = yb(4);
            dx1dswitch_b = yb(5);
            dx2dswitch_b = yb(6);
            dp1dswitch_b = yb(7);
            dp2dswitch_b = yb(8);

            res = [x1_a;
                   x2_a - 2
                   dx1dswitch_a;
                   dx2dswitch_a;
                   p1_b - x1_b + 4;
                   p2_b - x2_b + 2;
                   dp1dswitch_b - dx1dswitch_b;
                   dp2dswitch_b - dx2dswitch_b];
        end
    end
        
    function u = evalu(x, y)
        % Gets u at the sampled \tau 
        u = zeros(1, size(x, 2));
        for i = 1:size(x, 2)
            tau = x(1, i);
            
            p1 = y(3, i);
            p2 = y(4, i);
            
            if 0 <= tau && tau < 1
                % \tau \in [0, 1)
                u(1, i) = -p1 - p2;
            else
                % \tau \in [1, 2]
                u(1, i) = p2 - 2 * p1;
            end
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
                t(1, i) = tswitch + (2 - tswitch) * (tau - 1);
            end
        end
    end

    function J = evalJ(x, y, tswitch)
        t = evalt(x, tswitch);
        u = evalu(x, y);
        
        % Add the final cost.
        J = 0.5 * (y(1, end) - 4)^2 + 0.5 * (y(2, end) - 2)^2;
        
        % Add the running cost.
        for i = 1:size(t, 2) - 1
            dt = t(1, i + 1) - t(1, i);
            J = J + 0.5 * ((y(2, i) - 2)^2 + u(1, i)^2) * dt;
        end
    end

    function dJ = evaldJ(x, y, tswitch)
        % Add the contribution from the final cost.
        dJ = (y(1, end) - 4) * y(5, end) + (y(2, end) - 2) * y(6, end);
        
        % Add the contribution from the running cost.
        for i = 1:size(x, 2) - 1
            tau = x(1, i);
            dtau = x(1, i + 1) - x(1, i);
            
            x1 = y(1, i);
            x2 = y(2, i);
            p1 = y(3, i);
            p2 = y(4, i);
            dx1dswitch = y(5, i);
            dx2dswitch = y(6, i);
            dp1dswitch = y(7, i);
            dp2dswitch = y(8, i);
            
            if 0 <= tau && tau < 1
                % \tau \in [0, 1)
                u = -p1 - p2;
                dJ = dJ + dtau * (0.5 * ((x2 - 2)^2 + u^2) + tswitch * ((x2 - 2) * dx2dswitch + u * (-dp1dswitch - dp2dswitch)));
            else
                % \tau \in [1, 2]
                u = p2 - 2 * p1;
                dJ = dJ + dtau * (-0.5 * ((x2 - 2)^2 + u^2) + (2 - tswitch) * ((x2 - 2) * dx2dswitch + u * (dp2dswitch - 2 * dp1dswitch)));
            end
        end
    end
end