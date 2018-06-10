function test_cartpole_dynamics()

    % Draw the cartpole.
    cart_width = 0.6;
    cart_height = 0.3;
    cart_center = [0; 0.25];
    pole_length = 0.5;
    cart_mass = 10;
    pole_mass = 1;
    
    q = [0; pi - 0.15];
    dq = [0; 0];
    u = [0; 0];
    d = 0;
    figure(1);
    
    axis([-2.5, 2.5, -1, 2.5]);
    
    lc1 = line([q(1) + cart_center(1) - 0.5 * cart_width, q(1) + cart_center(1) - 0.5 * cart_width], [cart_center(2) - 0.5 * cart_height, cart_center(2) + 0.5 * cart_height]);
    lc2 = line([q(1) + cart_center(1) - 0.5 * cart_width, q(1) + cart_center(1) + 0.5 * cart_width], [cart_center(2) + 0.5 * cart_height, cart_center(2) + 0.5 * cart_height]);
    lc3 = line([q(1) + cart_center(1) + 0.5 * cart_width, q(1) + cart_center(1) + 0.5 * cart_width], [cart_center(2) + 0.5 * cart_height, cart_center(2) - 0.5 * cart_height]);
    lc4 = line([q(1) + cart_center(1) + 0.5 * cart_width, q(1) + cart_center(1) - 0.5 * cart_width], [cart_center(2) - 0.5 * cart_height, cart_center(2) - 0.5 * cart_height]);
    
    lp1 = line([q(1) + cart_center(1), q(1) + cart_center(1) + pole_length * sin(q(2))], [cart_center(2) + 0.5 * cart_height, cart_center(2) + 0.5 * cart_height + pole_length * cos(q(2))]);          
    
%     forward_dynamics(q, dq, ddq, u, cart_mass, pole_mass, pole_length)
    
    dt = 0.01;
    T = 1000;
    for i = 1:T
        % Update the coordinates of the cartpole.
        ddq = forward_dynamics(q, dq, u, cart_mass, pole_mass, pole_length);
        dq = dq + dt * ddq;
        q = q + dt * dq;
        
        % Update the visualization of the cartpole.
        set(lc1, 'Xdata', [q(1) + cart_center(1) - 0.5 * cart_width, q(1) + cart_center(1) - 0.5 * cart_width]);
        set(lc2, 'Xdata', [q(1) + cart_center(1) - 0.5 * cart_width, q(1) + cart_center(1) + 0.5 * cart_width]);
        set(lc3, 'Xdata', [q(1) + cart_center(1) + 0.5 * cart_width, q(1) + cart_center(1) + 0.5 * cart_width]);
        set(lc4, 'Xdata', [q(1) + cart_center(1) + 0.5 * cart_width, q(1) + cart_center(1) - 0.5 * cart_width]);
        set(lp1, 'Xdata', [q(1) + cart_center(1), q(1) + cart_center(1) + pole_length * sin(q(2))]);
        set(lp1, 'Ydata', [cart_center(2) + 0.5 * cart_height, cart_center(2) + 0.5 * cart_height + pole_length * cos(q(2))]);
        pause(dt);
    end

    function ddq = forward_dynamics(q, dq, u, mc, mp, l)
        g = 9.81;
        
%         ddq = [(1/(mc + mp)) * ((d + u(1)) - mp * l * (ddq(2) * cos(q(2)) - dq(2)^2 * sin(q(2))));
%                u(2)/(l * mp) + (g/l) * sin(q(2)) - (ddq(1)/l) * cos(q(2))];
        ddq = [(l * mp * sin(q(2))*dq(2)^2 + d + u(1) - u(2) * cos(q(2)) - g * mp * cos(q(2)) * sin(q(2)))/(mc + mp - mp*cos(q(2))^2);
               (mc * u(2) + mp * u(2) - d * mp * cos(q(2)) - mp * u(1) * cos(q(2)) + g * mp^2 * sin(q(2)) + g * mc * mp * sin(q(2)) - dq(2)^2 * l * mp^2 * cos(q(2)) * sin(q(2)))/(l * mp * (mc + mp - mp * cos(q(2))^2))];
    end
end