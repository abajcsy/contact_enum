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
        for i in range(len(u)):
            traj_opt.AddConstraintToAllKnotPoints(limits_low[i] <= u[i])
            traj_opt.AddConstraintToAllKnotPoints(u[i] <= limits_upp[i])

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

    def plan(self):
        print("not implemented!")
        # need to implement loop over mode sequences
        return None

    def _plan_hybrid_traj(self):
        # need to implement
        return None

