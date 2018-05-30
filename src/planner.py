#/usr/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt
import itertools

from pydrake.all import (
    DirectCollocation,
    PiecewisePolynomial,
    SolutionResult
)

class Planner:
    def __init__(self, plant, context, running_cost, signed_dist_funcs=[]):
        self.plant = plant
        self.context = context
        self.running_cost = running_cost
        self.signed_dist_funcs = signed_dist_funcs
        self.opt_params = {'num_time_samples': 21,
                           'minimum_timestep': 0.1,
                           'maximum_timestep': 0.4}

    def _solve_traj_opt(self, initial_state, final_state=None, duration_bounds=None, verbose=False):
        '''Finds a trajectory from an initial state, optionally to a final state.
        
        Args: 
            initial_state (tuple): the initial state
            final_state (tuple): the final state (default to None, final state unconstrained)
            duration (tuple): the min and max duration of the trajectory (default to None, 
                              no duration constraints)
            verbose (bool):

        Returns:
            pydrake.trajectories.PiecewisePolynomial: the planned trajectory
            float: the cost of the planned trajectory

        Raises:
            RuntimeError: raised if the optimization fails
        '''

        print("Initial state: {}\nFinal state: {}\nMin duration: {} s\nMax duration: {} s".format(
            initial_state, final_state, duration_bounds[0], duration_bounds[1]))

        traj_opt = DirectCollocation(self.plant, self.context, 
                                     self.opt_params['num_time_samples'],
                                     self.opt_params['minimum_timestep'],
                                     self.opt_params['maximum_timestep'])

        traj_opt.AddEqualTimeIntervalsConstraints()

        # Add bounds on the total duration of the trajectory
        if duration_bounds:
            traj_opt.AddDurationBounds(duration_bounds[0], duration_bounds[1])

        # TODO make input limits a paramter
        limits_low = [-15., -15.]
        limits_upp = [15., 15.]

        x = traj_opt.state()
        u = traj_opt.input()
        # for i in range(len(u)):
        #     traj_opt.AddConstraintToAllKnotPoints(limits_low[i] <= u[i])
        #     traj_opt.AddConstraintToAllKnotPoints(u[i] <= limits_upp[i])

        for signed_dist_func in self.signed_dist_funcs:
            traj_opt.AddConstraintToAllKnotPoints(signed_dist_func(x) >= 0)

        traj_opt.AddRunningCost(self.running_cost(x, u))

        # Add initial and final state constraints
        traj_opt.AddBoundingBoxConstraint(initial_state, initial_state,
                                          traj_opt.initial_state())

        if final_state:
            traj_opt.AddBoundingBoxConstraint(final_state, final_state,
                                              traj_opt.final_state())

            initial_x_trajectory = PiecewisePolynomial.FirstOrderHold([0., 0.4 * 21],
                                                                      np.column_stack((initial_state,
                                                                                       final_state)))
            traj_opt.SetInitialTrajectory(PiecewisePolynomial(), initial_x_trajectory)
        else:
            initial_x_trajectory = PiecewisePolynomial.FirstOrderHold([0., 0.4 * 21],
                                                                      np.column_stack((initial_state,
                                                                                       initial_state)))
            traj_opt.SetInitialTrajectory(PiecewisePolynomial(), initial_x_trajectory)

        result = traj_opt.Solve()

        if result != SolutionResult.kSolutionFound:
            raise RuntimeError('Direct collocation failed from initial state {}!'.format(initial_state))

        state_samples = traj_opt.GetStateSamples()
        input_samples = traj_opt.GetInputSamples()
        total_cost = 0.
        for k in range(state_samples.shape[1]):
            total_cost += self.running_cost(state_samples[:, k], input_samples[:, k])
            if verbose:
                for i, phi in enumerate(self.signed_dist_funcs):
                    print("\tsigned dist {}: {}".format(i, signed_dist_func(state_samples[:, k])))

        return traj_opt.ReconstructStateTrajectory(), total_cost

    def _enumerate_modes(self, tmin, T, dt, dc):
        """
        Returns a set of contact points and contact times. The contact points
        are defined by the y-coordinate of contact along the wall (since the 
        x-coordinate is fixed). 'None' for contact or time represents a 
        non-contact mode.
        """
        # TODO Need to make this general for arbitrary contact surfaces
        #cmin = 0.0
        cmin = 0.4
        cmax = 0.8
        #cmax = 0.75
        contacts = np.linspace(cmin,cmax,cmax/dc)
        times = np.linspace(tmin,T,T/dt)
        
        modes = list(itertools.product(contacts, times)) + [(None, None)]
        return modes

    def _plan_hybrid_traj(self, initial_state, c, t, T):
        """
        Returns a trajectory and associated cost (x_traj, total_cost) for a 
        given hybrid mode. If c and t are None, then we are simply planning a 
        contact-free trajectory. 
        """
        if c is None:
            # TODO raise exception on failure
            traj_x, total_cost = self._solve_traj_opt(initial_state, None, (T, T))

            return traj_x, None, total_cost

        # in contact mode, plan with contact-constrained final state
        x_wall = 2.0
        pole_len = 0.5
        cart_height = 0.4
        theta = math.acos((c - cart_height) / pole_len)
        x = x_wall - math.sqrt(pole_len**2 - (c - cart_height)**2)

        # TODO Need to address the final velocities = 0?
        final_state = (round(x, 8), round(theta, 8), 0., 0.) # TODO why is there a rounding error?

        # TODO raise exception on failure
        x_traj_nc, cost_nc = self._solve_traj_opt(initial_state, final_state, (t, t))

        # append zero control and x_traj_t..T constant at final_state
        if t >= T:
            x_traj_c = None
        else:
            x_traj_c = PiecewisePolynomial.FirstOrderHold([t, T],
                                                          np.column_stack((final_state,
                                                                           final_state)))
        cost_c = 0. # TODO compute this cost 

        return x_traj_nc, x_traj_c, cost_nc + cost_c

    def plan(self, initial_state, tmin, T, dt, dc):
        """
        Returns the minimum cost hybrid trajectory (either in contact
        or not in contact).
        tmin >= 0 is when to switch into contact mode, T is final time
        """
        #(contacts, times) = self._enumerate_modes(tmin, T, dt, dc)
        modes = self._enumerate_modes(tmin, T, dt, dc)
        #num_trajs = len(contacts)*len(times)
        num_trajs = len(modes)
        trajs = [None]*num_trajs
        costs = np.zeros(num_trajs)
        idx = 0
        # for c in contacts:
        #     for t in times:
        #         x_traj_nc, x_traj_c, total_cost = self._plan_hybrid_traj(initial_state, c, t, T)
        #         trajs[idx] = (x_traj_nc, x_traj_c)
        #         costs[idx] = total_cost
        #         idx += 1
        
        for (c, t) in modes:
            x_traj_nc, x_traj_c, total_cost = self._plan_hybrid_traj(initial_state, c, t, T)
            trajs[idx] = (x_traj_nc, x_traj_c)
            costs[idx] = total_cost
            idx += 1
            
        min_idx = np.argmin(costs)
        return trajs[min_idx]

