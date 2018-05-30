#/usr/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt

from pydrake.multibody.rigid_body_tree import (
    RigidBodyTree,
    FloatingBaseType,
    AddModelInstanceFromUrdfFile
)
from pydrake.multibody.rigid_body_plant import (
    RigidBodyPlant
)
from underactuated import (
    PlanarRigidBodyVisualizer
)

from planner import Planner


def main():
    tree = RigidBodyTree()
    AddModelInstanceFromUrdfFile("resources/cartpole2.urdf", 
                                 FloatingBaseType.kFixed,
                                 None,
                                 tree)
    AddModelInstanceFromUrdfFile("resources/wall.urdf", 
                                 FloatingBaseType.kFixed,
                                 None,
                                 tree)

    plant = RigidBodyPlant(tree)
    context = plant.CreateDefaultContext()

    #def running_cost(x, u):
    #    return 10. * (u[0]**2 + u[1]**2)
    
    input_weight = 1.
    goal_weight = 30.
    def running_cost(x, u):
        return input_weight * (u[0]**2 + u[1]**2) + goal_weight * (x[0] - 2.)**2

    import pydrake.symbolic as sym

    def signed_dist_func(x):
        pole_x = 0.5 * sym.cos(x[1]) + x[0]
        #pole_x = 0
        return 2. - pole_x

    planner = Planner(plant, context, running_cost, signed_dist_func)

    #initial_state = (0., math.pi, 0., 0.)
    initial_state = (0., 0., 0., 0.)
    #final_state = (1., 3*math.pi/2., 0., 0.)
    final_state = None
    x_traj, total_cost = planner._solve_traj_opt(initial_state, final_state)

    vis = PlanarRigidBodyVisualizer(tree, xlim=[-2.5, 2.5], ylim=[-1, 2.5])
    ani = vis.animate(x_traj, repeat=True)

    plt.show()    


if __name__ == '__main__':
    main()
