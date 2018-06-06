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
    def __init__(self, plant, context, running_cost, final_cost, signed_dist_funcs=[], final_state_constraint=None):
        self.plant = plant
        self.context = context
        self.running_cost = running_cost
        self.final_cost = final_cost
        self.signed_dist_funcs = signed_dist_funcs
        self.final_state_constraint = final_state_constraint
        self.opt_params = {'num_time_samples': 21,
                           'minimum_timestep': 0.1,
                           'maximum_timestep': 0.4}

    def _solve_traj_opt(self, initial_state, constrain_final_state=True, duration_bounds=None, d=0.0, verbose=False):
        '''Finds a trajectory from an initial state, optionally to a final state.
        
        Args: 
            initial_state (tuple): the initial state
            final_state (tuple): the final state (default to None, final state unconstrained)
            duration (tuple): the min and max duration of the trajectory (default to None, 
                              no duration constraints)
            d (float): constant disturbance force
            verbose (bool): enables/disables verbose output

        Returns:
            pydrake.trajectories.PiecewisePolynomial: the planned trajectory
            float: the cost of the planned trajectory

        Raises:
            RuntimeError: raised if the optimization fails
        '''

        print("Initial state: {}\nFinal state: {}\nMin duration: {} s\nMax duration: {} s".format(
            initial_state, constrain_final_state, duration_bounds[0], duration_bounds[1]))

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
        t = traj_opt.time()
        # TODO assuming disturbance is at the last index
        for i in range(len(u) - 1):
            traj_opt.AddConstraintToAllKnotPoints(limits_low[i] <= u[i])
            traj_opt.AddConstraintToAllKnotPoints(u[i] <= limits_upp[i])

        traj_opt.AddConstraintToAllKnotPoints(u[len(u) - 1] == d)

        for signed_dist_func in self.signed_dist_funcs:
            traj_opt.AddConstraintToAllKnotPoints(signed_dist_func(x) >= 0)

        traj_opt.AddRunningCost(traj_opt.timestep(0) * self.running_cost(x, u, t))

        traj_opt.AddFinalCost(self.final_cost(x, u, t))

        # Add initial and final state constraints
        traj_opt.AddBoundingBoxConstraint(initial_state, initial_state,
                                          traj_opt.initial_state())

        if self.final_state_constraint and constrain_final_state:
            traj_opt.AddConstraint(self.final_state_constraint(traj_opt.final_state()) == 0)

        # # TODO this is redundant with the final state equality constraint above
        # if final_state:
        #     traj_opt.AddBoundingBoxConstraint(final_state, final_state,
        #                                       traj_opt.final_state())

        #     initial_x_trajectory = PiecewisePolynomial.FirstOrderHold([0., 0.4 * 21],
        #                                                               np.column_stack((initial_state,
        #                                                                                final_state)))
        #     traj_opt.SetInitialTrajectory(PiecewisePolynomial(), initial_x_trajectory)
        # else:
        #     initial_x_trajectory = PiecewisePolynomial.FirstOrderHold([0., 0.4 * 21],
        #                                                               np.column_stack((initial_state,
        #                                                                                initial_state)))
        #     traj_opt.SetInitialTrajectory(PiecewisePolynomial(), initial_x_trajectory)

        initial_x_trajectory = PiecewisePolynomial.FirstOrderHold([0., 0.4 * 21],
                                                                  np.column_stack((initial_state,
                                                                                   initial_state)))
        traj_opt.SetInitialTrajectory(PiecewisePolynomial(), initial_x_trajectory)

        result = traj_opt.Solve()

        if result != SolutionResult.kSolutionFound:
            raise RuntimeError('Direct collocation failed from initial state {}!'.format(initial_state))

        state_samples = traj_opt.GetStateSamples()
        input_samples = traj_opt.GetInputSamples()
        time_samples = traj_opt.GetSampleTimes()

        # for debugging
        hs = [time_samples[i+1] - time_samples[i] for i in range(len(time_samples)) if i < len(time_samples) - 1]
        #print(hs)

        total_cost = 0.
        for k in range(state_samples.shape[1]):
            total_cost += (hs[0] * 
                           self.running_cost(state_samples[:, k], 
                                             input_samples[:, k], 
                                             time_samples[k]))
            if verbose:
                for i, phi in enumerate(self.signed_dist_funcs):
                    print("\tsigned dist {}: {}".format(i, signed_dist_func(state_samples[:, k])))

        if verbose:
            print("Total cost is {}".format(total_cost))

            u_traj = traj_opt.ReconstructInputTrajectory()
            times = np.linspace(u_traj.start_time(), u_traj.end_time(), 100)
            u_lookup = np.vectorize(lambda t: u_traj.value(t)[0])
            u_values = u_lookup(times)

            plt.figure()
            plt.plot(times, u_values)
            plt.xlabel('time (seconds)')
            plt.ylabel('force (Newtons)')

            plt.show()

        return traj_opt.ReconstructStateTrajectory(), total_cost

    def _enumerate_modes(self, tmin, T, dt, dc):
        """
        Returns a set of contact points and contact times. The contact points
        are defined by the y-coordinate of contact along the wall (since the 
        x-coordinate is fixed). 'None' for contact or time represents a 
        non-contact mode.
        """
        modes = list(np.linspace(tmin, T, T/dt)) + [None]
        return modes

    def _plan_hybrid_traj(self, initial_state, t, T, d):
        """
        Returns a trajectory and associated cost (x_traj, total_cost) for a 
        given hybrid mode. If c and t are None, then we are simply planning a 
        contact-free trajectory. 
        """
        try:
            if t is None:
                print("Planning hybrid traj without final state constraint")
                traj_x, total_cost = self._solve_traj_opt(initial_state, False, (T, T), d)

                return traj_x, None, total_cost

            x_traj_nc, cost_nc = self._solve_traj_opt(initial_state, True, (t, t), d)

            # append zero control and x_traj_t..T constant at final_state
            if t >= T:
                x_traj_c = None
            else:
                # TODO fix this
                x_traj_c = PiecewisePolynomial.FirstOrderHold([t, T],
                                                              np.column_stack((initial_state,
                                                                               initial_state)))
            cost_c = 0. # TODO compute this cost 

            return x_traj_nc, x_traj_c, cost_nc + cost_c
        except RuntimeError as e:
            raise e

    def plan(self, initial_state, tmin, T, dt, dc, d):
        """
        Returns the minimum cost hybrid trajectory (either in contact
        or not in contact).
        tmin >= 0 is when to switch into contact mode, T is final time
        """
        modes = self._enumerate_modes(tmin, T, dt, dc)
        num_trajs = len(modes)
        trajs = [(None, None)] * num_trajs
        costs = np.zeros(num_trajs)
        idx = 0
        
        #for k, (c, t) in enumerate(modes):
        for k, t in enumerate(modes):
            try:
                x_traj_nc, x_traj_c, total_cost = self._plan_hybrid_traj(initial_state, t, T, d)
                trajs[k] = (x_traj_nc, x_traj_c)
                costs[k] = total_cost
                print("Cost of trajectory {} of {} is {}".format(k, len(modes) - 1, costs[k]))
            except RuntimeError as e:
                print("Failed to plan a hybrid trajectory!")
                trajs[k] = (None, None)
                costs[k] = 1e9
            
        min_idx = np.argmin(costs)

        in_contact = trajs[min_idx][1] is not None
        print("Choosing trajectory at index {} with cost {} (in contact? {})".format(min_idx, costs[min_idx], in_contact))

        return trajs[min_idx]

