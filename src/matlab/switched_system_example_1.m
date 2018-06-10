classdef switched_system_example_1
    %SWITCHED_SYSTEM_EXAMPLE_1 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t0
        tf
    end
    
    methods
        function obj = switched_system_example_1(t0, tf)
            obj.t0 = t0;
            obj.tf = tf;
        end
        
        %% Returns dynamics of switched system:   
        %   dx(t)/dt = f_i(x(z), u(z)), t_{i-1} <= t < t_i
        function [dx, xnext] = dyn(obj, x, u, t, dt, stimes)
            if t >= 0 && t < stimes(1)
                dx = [x(1) + u*sin(x(1)); -x(2) + u*cos(x(2))];
            elseif t >= stimes(1) && t < stimes(2)
                dx = [x(2) + u*sin(x(2)); -x(1) - u*cos(x(1))];
            else
                dx = [-x(1) - u*sin(x(1)); x(2) + u*cos(x(2))];
            end
            % do first order Euler integration
            xnext = x' + dt*dx;
        end
        
        %% Returns time-normalized dynamics of switched system
        %   dx(z)/dz = (t_i - t_i-1)f_i(x(z), u(z)), i-1 <= z < i
        %   t = (t_i - t_i-1)*(z-i) + t_i
        function [dx, xnext, t] = dyn_z(obj, x, u, z, dz, stimes)
            % convert time into normalized time variable, z
            if z >= 0 && z < 1
                dx = (stimes(1) - obj.t0)*[x(1) + u*sin(x(1)); -x(2) + u*cos(x(2))];
                t = (stimes(1) - obj.t0)*(z-1)+stimes(1);
            elseif z >= 1 && z < 2
                dx = (stimes(2) - stimes(1))*[x(2) + u*sin(x(2)); -x(1) - u*cos(x(1))];
                t = (stimes(2) - stimes(1))*(z-2)+stimes(2);
            else
                dx = (obj.tf - stimes(2))*[-x(1) - u*sin(x(1)); x(2) + u*cos(x(2))];
                t = (obj.tf - stimes(2))*(z-3)+obj.tf;
            end
            % TODO: what should dz be?
            xnext = x' + dz*dx;
        end

%         %% Returns linearized dynamics about (xbar, ubar)
%         %   d(deltax)/dt = A_i(t)*deltax + B_i(t)*deltau
%         function [dx, xnext] = dyn_lin(obj, x, xbar, u, ubar, t, dt, stimes)
%             [A, B] = linearize(xbar, ubar, t, stimes);
%             deltax = x-xbar;
%             deltau = u-ubar;
%             dx = A*deltax + B*deltau;
%             xnext = x' + dt*dx;
%         end
        
        %% Returns linearized and normalized dynamics about (xbar, ubar)
        %   d(deltax)/dz = (t_i - t_{i-1})*(A_i(z)*deltax + B_i(z)*deltau)
        function [dx, xnext] = dyn_z_lin(obj, x, xbar, u, ubar, z, dz, stimes)
            [A, B] = linearize(xbar, ubar, z);
            deltax = x-xbar;
            deltau = u-ubar;
            dx = A*deltax + B*deltau;
            if z >= 0 && z < 1
                dx = (stimes(1) - obj.t0)*dx;
            elseif z >= 1 && z < 2
                dx = (stimes(2) - stimes(1))*dx;
            else
                dx = (obj.tf - stimes(2))*dx;
            end
            xnext = x' + dz*dx;
        end
        
        %% Linearizes the switched system about nominal (xbar, ubar) 
        function [A, B] = linearize(obj, xbar, ubar, z)
            if z >= 0 && z < 1
                A = [1 + ubar*cos(xbar(1)), 0; 
                    -1 - ubar*sin(xbar(2)), 0];
                B = [sin(xbar(1)); cos(xbar(2))];
            elseif z >= 1 && z < 2
                A = [0, 1 + ubar*cos(xbar(2)); 
                    -1 + ubar*sin(xbar(1)), 0];
                B = [sin(xbar(2)); -cos(xbar(1))];
            else
                A = [-1 - ubar*cos(xbar(1)), 0; 
                    0, 1 - ubar*sin(xbar(2))];
                B = [-sin(xbar(1)); cos(xbar(2))];
            end
        end
    end
end

