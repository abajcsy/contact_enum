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
        end
        
        %% Solves the SLQ problem given initial state and control rollout
        function [x, u, dS, dsvec, dsscalar] = run(obj, x0, uinit, stimes)
            ubar = uinit;
            obj.stimes = stimes;
            num_pts = (length(stimes)+1)/obj.dz; 
            
            iter = 1;
            while iter <= obj.max_iter % or norm(lvec) < lmin
                % forward integrate the system dynamics using x0 and uinit
                xbar = obj.rollout(x0, ubar);
                if iter == 1 
                    x = xbar; 
                    u = ubar;
                end
                dx = (x - xbar); % note: for first rollout, dx = 0
                
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
                
                W = cell(1,num_pts);
                wvec = cell(1,num_pts);
                wscalar = cell(1,num_pts);
                
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
                    
                    W{zidx+1} = Q + Amats{zidx+1}'*S{zidx+1} + S{zidx+1}*Amats{zidx+1} - L{zidx+1}'*R*L{zidx+1};
                    wvec{zidx+1} = qvec + Amats{zidx+1}'*svec{zidx+1} - L{zidx+1}'*R*lvec{zidx+1};
                    wscalar{zidx+1} = qscalar - 0.5*alpha*(2-alpha)*lvec{zidx+1}'*R*lvec{zidx+1};

                    % get the real time interval for current z value
                    z = zvals{zidx};
                    [~, ~, t_i, t_i1] = obj.convert_z(z);
                    
                    % integrate S,svec,sscalar  
                    %   -dS/dz = (t_i-t_{i-1})*W, 
                    %   -dsvec/dz = (t_i-t_{i-1})*wvec, 
                    %   -dsscalar/dz = (t_i-t_{i-1})*wscalar
                    S{zidx} = S{zidx+1} - (t_i-t_i1)*W{zidx+1};
                    svec{zidx} = svec{zidx+1} - (t_i-t_i1)*wvec{zidx+1};
                    sscalar{zidx} = sscalar{zidx+1} - (t_i-t_i1)*wscalar{zidx+1};
                    
                    % TODO: armijo line search for the optimal alpha with policy
                                    
                    % update control law for current time
                    ubar(zidx+1) = ubar(zidx+1) + alpha*lvec{zidx+1}' + L{zidx+1}*dx(zidx+1,:)';
                end
                
                % compute final L, lvec
                [~,~,~,~, rvec, R, P]= ...
                        obj.lin_cost(x(zidx), xbar(zidx), u(zidx), ubar(zidx));
                L{zidx} = -inv(R)*(P' + Bmats{zidx}'*S{zidx});
                lvec{zidx} = -inv(R)*(rvec + Bmats{zidx}'*svec{zidx});

                W{zidx} = Q + Amats{zidx+1}'*S{zidx+1} + S{zidx+1}*Amats{zidx+1} - L{zidx+1}'*R*L{zidx+1};
                wvec{zidx} = qvec + Amats{zidx+1}'*svec{zidx+1} - L{zidx+1}'*R*lvec{zidx+1};
                wscalar{zidx} = qscalar - 0.5*alpha*(2-alpha)*lvec{zidx+1}'*R*lvec{zidx+1};
                
                % update final control signal
                ubar(zidx) = ubar(zidx) + alpha*lvec{zidx}' + L{zidx}*dx(zidx,:)';
                
                % store the prior rollout as the next candidate state seq
                x = xbar;
                u = ubar;
                
                % update iterate
                iter = iter +1;
            end
            
            % compute the partials wrt switching time
            %   cell array rows correspond to switching time
            %   cell array columns correspond to partials at time z
            dS = cell(2,num_pts);
            dsvec = cell(2,num_pts);
            dsscalar = cell(2,num_pts);
            
            % set the final conditions for dS,dsvec,dsscalar
            dS{1,num_pts} = zeros(2);
            dS{2,num_pts} = zeros(2);
            dsvec{1,num_pts} = Qf*partialx(num_pts);
            dsvec{2,num_pts} = Qf*partialx(num_pts);
            dsscalar{1,num_pts} = qfvec'*partialx(num_pts);
            dsscalar{2,num_pts} = qfvec'*partialx(num_pts);
            
            dR = 0;
            drvec = cell(1,num_pts);
            dP = [0;0];
            
            dQ = zeros(2); %TODO: THIS IS PROBS WRONG DIM
            dqvec = cell(1,num_pts);
            dqscalar = cell(1,num_pts);
            
            for j=1:2 % which switching time we are taking partials w.r.t.
                % store partial of x and u w.r.t. switching time
                partialx = cell(1,num_pts);
                partialu = cell(1,num_pts);
                
                partialx{1} = [0;0];
                partialu{1} = -L{1}*partialx{1} + dL{1}*(x-xbar)+dlvec{1};

                % store the partials of L and lvec 
                dL = cell(1,num_pts);
                dlvec = cell(1,num_pts);
                for zidx=(num_pts-1):-1:1
                    % get the real time interval for current z value
                    z = zvals{zidx};
                    [~, i, t_i, t_i1] = obj.convert_z(z);
                    
                    %delta_{i,j} = 1 if i=j, and 0 otherwise
                    delta_ij = i==j;
                    delta_i1j = (i-1)==j;

                    drvec{zidx+1} = P'*partialx{zidx+1} + R*partialu{zidx+1};
                    dqscalar{zidx+1} = qvec'*partialx{zidx+1} + r'*partialu{zidx+1};
                    dqvec{zidx+1} = Q*partialx{zidx+1} + P*partialu{zidx+1};
                    
                    % compute partial of switch time for L and lvec
                    dL{zidx+1} = -inv(R)*Bmats{zidx+1}'*dS{j,zidx+1};
                    dlvec{zvec+1} = -inv(R)*(drvec{zidx+1}+Bmats{zidx+1}'*dsvec{j,zidx+1});
          
                    % compute partial of u
                    partialu{zidx+1} = -L{zidx+1}*partialx{zidx+1}*dL{zidx+1} + dlvec{zidx+1};
                    
                    % compute partial of x
                    partialx{zidx} = partialx{zidx} ...
                        - obj.dz*(deltaij - deltai1j)*obj.plant.dyn_z(x(zidx), u(zidx), z, obj.dz, obj.stimes) ...
                        + (t_i - t_i1)*(Amats{zidx+1}*partialx{zidx+1} + Bmats{zidx+1}*partialu{zidx+1});
                    
                    dW = Amats{zidx+1}'*dS{j,zidx+1} + dS{j,zidx+1}*A{zidx+1} - dL'*R*L{zidx+1} - L{zidx+1}'*R*dL;
                    dwvec = dqvec + A{zidx+1}'*dsvec{j,zidz+1} - dL'*R*lvec - L{zidx+1}'*R*dlvec;
                    dwscalar = dqscalar-0.5*alpha*(2-alpha)*(dlvec'*R*lvec)-lvec'*R*dlvec;

                    %-dS = (delta_{i,j} - delta_{i-1,j})*W+(ti-ti-1)*dW
                    %-dsvec = (delta_{i,j} - delta_{i-1,j})*wvec+(ti-ti-1)*dwvec
                    %-dsscalar = (delta_{i,j} - delta_{i-1,j})*wscalar+(ti-ti-1)*dwscalar
                    dS{j,zidx} = dS{zidx+1} - ((delta_ij - delta_i1j)*W{zidx+1} + (t_i-t_i1)*dW{zidx+1});
                    dsvec{j,zidx} = dsvec{zidx+1} - ((delta_ij - delta_i1j)*wvec{zidx+1} + (t_i-t_i1)*dwvec{zidx+1});
                    dsscalar{j,zidx} = dsscalar{zidx+1} - ((delta_ij - delta_i1j)*wscalar{zidx+1} + (t_i-t_i1)*dwscalar{zidx+1});
                end
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
        function [t, i, t_i, t_i1] = convert_z(obj, z)
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

