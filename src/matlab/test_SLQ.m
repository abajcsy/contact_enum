% these are optimal switching times base on (Xu, 2004)!
stimes = [0.2262; 1.0176];
xg = [1; -1];
t0 = 0;
tf = 3;
dt = 0.01;
dz = 0.01;

% create the plant and the controller
system = switched_system_example_1(t0, tf);
controller = SLQ(system, dt, dz, t0, tf, xg);

% define initial condition and initial control sequence
x0 = [2; 3];
uinit = 0*ones((3-0)/dz, 1);

% solve the SLQ problem
[ubar, dS, dsvec, dsscalar] = controller.run(x0, uinit, stimes);

% get gradJ w.r.t switching time vector
gradJ = [dsscalar{1,1}; dsscalar{2,1}];

% use gradient method to update switching times
alpha = 0.5;
stimes = stimes + alpha*gradJ;
