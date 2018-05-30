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

        return traj_opt.ReconstructStateTrajectory()

    def plan(self):
        print("not implemented!")
        return None

