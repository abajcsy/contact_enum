#/usr/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt

from pydrake.all import (
    DirectCollocation,
    PiecewisePolynomial,
    SolutionResult
)

class Planner:
    def __init__(self, plant, context, running_cost, signed_dist_func):
        self.plant = plant
        self.context = context
        self.running_cost = running_cost
        self.signed_dist_func = signed_dist_func
        self.opt_params = {'num_time_samples': 21,
                           'minimum_timestep': 0.1,
                           'maximum_timestep': 0.4}

    def _solve_traj_opt(self, initial_state, final_state=None):
        traj_opt = DirectCollocation(self.plant, self.context, 
                                     self.opt_params['num_time_samples'],
                                     self.opt_params['minimum_timestep'],
                                     self.opt_params['maximum_timestep'])

        traj_opt.AddEqualTimeIntervalsConstraints()

        # TODO make input limits a paramter
        limits_low = [-15., -15.]
        limits_upp = [15., 15.]

        x = traj_opt.state()
        u = traj_opt.input()
        for i in range(len(u)):
            traj_opt.AddConstraintToAllKnotPoints(limits_low[i] <= u[i])
            traj_opt.AddConstraintToAllKnotPoints(u[i] <= limits_upp[i])

        traj_opt.AddConstraintToAllKnotPoints(self.signed_dist_func(x) >= 0)

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

        result = traj_opt.Solve()
        assert(result == SolutionResult.kSolutionFound)

        state_samples = traj_opt.GetStateSamples()
        input_samples = traj_opt.GetInputSamples()
        total_cost = 0.
        for k in range(state_samples.shape[1]):
            total_cost += self.running_cost(state_samples[:, k], input_samples[:, k])
            print(self.signed_dist_func(state_samples[:, k]))

        print("Total cost: {}".format(total_cost))

        return traj_opt.ReconstructStateTrajectory(), total_cost


    def _enumerate_modes(self, tmin, T, dt, dc):
        """
        Returns a set of contact points and contact times. The contact points
        are defined by the y-coordinate of contact along the wall (since the 
        x-coordinate is fixed). 'None' for contact or time represents a 
        non-contact mode.
        """
        # TODO Need to make this general for arbitrary contact surfaces
        cmin = 0.0
        cmax = 0.75
        contacts = np.append(np.linspace(cmin,cmax,cmax/dc), np.array(None))
        times = np.append(np.linspace(tmin,T,T/dt), np.array(None))
        return (contacts, times)

    def _plan_hybrid_traj(self, initial_state, c, t):
        """
        Returns a trajectory and associated cost (x_traj, total_cost) for a 
        given hybrid mode. If c and t are None, then we are simply planning a 
        contact-free trajectory. 
        """
        final_state = None

        if c is not None:
            # in contact mode, plan with contact-constrained final state
            x_wall = 2.0
            pole_len = 0.5
            cart_height = 0.4
            theta = np.arcsin((c - cart_height)/pole_len)
            x = x_wall - np.sqrt(pole_len**2 - (c - cart_height)**2)
            # TODO Need to address the final velocities = 0?
            final_state = (x, theta, 0., 0.)

        x_traj, total_cost = planner._solve_traj_opt(initial_state, final_state)

        if c is not None:
            # append zero control and x_traj_t..T constant at final_state
            
        return x_traj, total_cost

    def plan(self, initial_state, tmin, T, dt, dc):
        """
        Returns the minimum cost hybrid trajectory (either in contact
        or not in contact).
        tmin >= 0 is when to switch into contact mode, T is final time
        """
        (contacts, times) = self._enumerate_modes(tmin, tmax, dt, dc)
        num_trajs = len(contacts)*len(times)
        trajs = [None]*num_trajs
        costs = np.zeros(num_trajs)
        idx = 0
        for c in contacts
            for t in times
                x_traj, total_cost = self._plan_hybrid_traj(initial_state, c, t)
                trajs[idx] = x_traj
                costs[idx] = total_cost
                idx += 1
        min_idx = np.argmin(costs)
        return trajs[min_idx]

