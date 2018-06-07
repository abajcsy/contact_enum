classdef SLQ
    %SLQ Implementation of Sequential Linear Quadratic control 
    %    based on the paper:
    %    Farshidian, F., Kamgarpour, M., Pardo, D., Buchli, J. (2017) 
    %       Sequential Linear Quadratic Optimal Control for Nonlinear 
    %       Switched Systems. International Federation of Automatic
    %       Control (IFAC). 
    
    properties
        plant
        dt
        dz
        t0
        tf
        xg
        num_states
        num_control
        convergence_thresh
        max_iter
        stimes
    end
    
    methods
        function obj = SLQ(plant, dt, dz, t0, tf, xg)
            %SLQ Construct an instance of this class
            %   Detailed explanation goes here
            obj.plant = plant;
            obj.dt = dt;
            obj.dz = dz;
            obj.xg = xg;
            obj.t0 = t0;
            obj.tf = tf;
            obj.num_states = length(xg);
            obj.num_control = 1;
            obj.max_iter = 1;
            obj.stimes = [];
            
            % final cost coefficients
            % Phi_f(x) = qfscalar + qfvec'*dx + 0.5*dx'*Qf*dx
            
            % running cost coefficients
            % L(z,x,u) = qscalar + dx'*qvec + du'*r + dx'*P*du + 0.5*dx'*Q*dx + 0.5*du'*R*du
        end
        
        %% Solves the SLQ problem given initial state and control rollout
        function [ubar, dS, dsvec, dsscalar] = run(obj, x0, uinit, stimes)
            ubar = uinit;
            obj.stimes = stimes;
            for ii=1:obj.max_iter
                % forward integrate the system dynamics using x0 and uinit
                xbar = obj.rollout(x0, ubar);
                if ii == 1 
                    x = xbar; 
                    u = ubar;
                end
                dx = (x - xbar); % note: for first rollout, dx = 0
                num_pts = obj.tf/obj.dt; % TODO: change this to dz?
                
                Amats = cell(1,num_pts);
                Bmats = cell(1,num_pts);
                zvals = cell(1,num_pts);
                z = 0;
                % compute the LQ approximation to the dynamics
                for idx=1:num_pts
                    zvals{idx} = z;
                    [A, B] = obj.plant.linearize(xbar(idx,:), ubar(idx,:), z);
                    Amats{idx} = A;
                    Bmats{idx} = B;
                    z = z + obj.dz;
                end
                
                alpha = 0.1; %TODO WHAT IF ARMIJO LINE SEARCH IS SENTIENT?
                
                % solve the Riccati equations to find the optimal control
                S = cell(1,num_pts);
                svec = cell(1,num_pts);
                sscalar = cell(1,num_pts);
                L = cell(1,num_pts);
                lvec = cell(1,num_pts);
                
                % compute the partials wrt switching time
                %   cell array rows correspond to switching time
                %   cell array columns correspond to partials at time z
                dS = cell(2,num_pts);
                dsvec = cell(2,num_pts);
                dsscalar = cell(2,num_pts);
                
                % get the quadraticized cost at final time
                [~, qfscalar, qfvec, Qf] = ...
                    obj.lin_cost_final(x(num_pts), xbar(num_pts));
                
                % set the final conditions for S,svec,sscalar, L,lvec
                S{num_pts} = Qf;
                svec{num_pts} = qfvec;
                sscalar{num_pts} = qfscalar;
                
                % solve the final value differential equations
                for zidx=(num_pts-1):-1:1
                    [~, qscalar, qvec, Q, rvec, R, P]= ...
                        obj.lin_cost(x(zidx), xbar(zidx), u(zidx), ubar(zidx));
                    
                    % compute L, lvec
                    L{zidx+1} = -inv(R)*(P' + Bmats{zidx+1}'*S{zidx+1});
                    lvec{zidx+1} = -inv(R)*(rvec + Bmats{zidx+1}'*svec{zidx+1});
                    
                    W = Q + Amats{zidx+1}'*S{zidx+1} + S{zidx+1}*Amats{zidx+1} - L{zidx+1}'*R*L{zidx+1};
                    wvec = qvec + Amats{zidx+1}'*svec{zidx+1} - L{zidx+1}'*R*lvec{zidx+1};
                    wscalar = qscalar - 0.5*alpha*(2-alpha)*lvec{zidx+1}'*R*lvec{zidx+1};

                    % get the real time interval for current z value
                    z = zvals{zidx};
                    [~, t_i, t_i1] = obj.convert_z(z);
                    
                    % integrate S,svec,sscalar  
                    %   -dS/dz = (t_i-t_{i-1})*W, 
                    %   -dsvec/dz = (t_i-t_{i-1})*wvec, 
                    %   -dsscalar/dz = (t_i-t_{i-1})*wscalar
                    S{zidx} = S{zidx+1} - (t_i-t_i1)*W;
                    svec{zidx} = svec{zidx+1} - (t_i-t_i1)*wvec;
                    sscalar{zidx} = sscalar{zidx+1} - (t_i-t_i1)*wscalar;
                end
                
                % compute final L, lvec
                [~,~,~,~, rvec, R, P]= ...
                        obj.lin_cost(x(zidx), xbar(zidx), u(zidx), ubar(zidx));
                L{zidx} = -inv(R)*(P' + Bmats{zidx+1}'*S{zidx+1});
                lvec{zidx} = -inv(R)*(rvec + Bmats{zidx+1}'*svec{zidx+1});
                
                % TODO: armijo line search for the optimal alpha with policy
                
                % update control law
                Lreshape = cellfun(@(x)reshape(x,[2,1]),L,'un',0);
                ubar(t) = ubar + alpha*cell2mat(lvec)' + cell2mat(Lreshape)*dx;
                
                % store the prior rollout as the next candidate state seq
                x = xbar;
                u = ubar;
                
                %V = sscalar + dx'*svec + 0.5*dx'*S*dx
            end
        end

        %% Forward integrate system dynamics using the latest estimation
        %   of the optimal control law with a fixed switching time vector.
        function x = rollout(obj, x0, u)
            tN = length(u);
            % TODO: should x be tN+1 size? (including x0)
            x = zeros(tN+1, obj.num_states); 
            x(1,:) = x0;
            z = 0;
            for idx=1:tN
                [~, xnext, ~] = obj.plant.dyn_z(x(idx,:), u(idx,:), z, obj.dz, obj.stimes);
                x(idx+1,:) = xnext;
                z = z + obj.dz;
            end
        end
        
        %% Compute total (quadraticized) cost for trajectory (x,u)
        %   J = Phi(x_f) + 
        %       \sum^3_i=1 \integral^i_{i-1} (t_i-t_{i-1})L(x,u) dz
        function cum_cost = cum_cost_z(obj, x, xbar, u, ubar)
            [cost_f,~,~,~] = obj.lin_cost_final(x(end), xbar(end));
            cum_cost = cost_f;
            idx = 1;
            for mode=1:3
                for z=mode-1:obj.dz:mode
                    [cost,~,~,~,~,~,~] = ...
                        obj.lin_cost(x(idx,:), xbar(idx,:), u(idx,:), ubar(idx,:));
                    cum_cost = cum_cost + cost;
                    idx = idx+1;
                end
            end
        end
        
        %% Quadraticized running cost function defined by
        %   L(x,u) = norm(u)^2
        function [cost, qscalar, qvec, Q, rvec, R, P]= lin_cost(obj, x, xbar, u, ubar)
            deltax = x-xbar;
            deltau = u-ubar;
            qscalar = ubar'*ubar;
            qvec = [0; 0];
            Q = [0, 0; 0, 0];
            rvec = 2*ubar;
            R = 2;
            P = [0; 0];
            cost = qscalar + deltax'*qvec + deltau'*rvec ...
                + deltax'*P*deltau + 0.5*deltax'*Q*deltax ...
                + 0.5*deltau'*R*deltau;
        end
        
        %% Quadraticized final cost function defined by
        %   Phi(x) = norm(x-xg)^2
        function [cost_f, qfscalar, qfvec, Qf] = lin_cost_final(obj, x, xbar)
            deltax = x-xbar;
            qfscalar = norm(xbar - obj.xg)^2;
            qfvec = 2*(xbar - obj.xg);
            Qf = 2*eye(2);
            cost_f = qfscalar + qfvec'*deltax + 0.5*deltax'*Qf*deltax;
        end
        
        %% Converts a z value to a (time, t_i, t_i-1)
        function [t, t_i, t_i1] = convert_z(obj, z)
            if z >= 0 && z < 1
               i = 1; 
               t_i = obj.stimes(1);
               t_i1 = obj.t0;
            elseif z >= 1 && z < 2
               i = 2;
               t_i = obj.stimes(2);
               t_i1 = obj.stimes(1);
            else
               i = 3;
               t_i = obj.tf;
               t_i1 = obj.stimes(2);
            end
            t = (t_i - t_i1)*(z-i)+t_i;
        end
    end
end

