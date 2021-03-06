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
from pydrake.systems.framework import (
    BasicVector
)
from underactuated import (
    PlanarRigidBodyVisualizer
)
from planner import Planner
import pydrake.symbolic as sym


def init_signed_dist_funcs():
    signed_dist_funcs = []

    wall_x = 2.
    pole_length = 0.5
    cart_width = 0.6
    cart_height = 0.4

    # Distance from the tip of the pole to the wall
    def dist_pole_wall(x):
        pole_x = pole_length * sym.sin(x[1]) + x[0]
        return wall_x - pole_x

    # Distance from the tip of the pole to the ground
    def dist_pole_ground(x):
        pole_y = pole_length * sym.cos(x[1]) + cart_height
        return pole_y

    # Distance from the front of the cart to the wall
    def dist_cart_wall(x):
        cart_x = x[0] + 0.5 * cart_width
        return wall_x - cart_x

    return [dist_pole_wall, dist_pole_ground, dist_cart_wall]


def test_planner():
    tree = RigidBodyTree()
    # AddModelInstanceFromUrdfFile("resources/cartpole2.urdf", 
    #                              FloatingBaseType.kFixed,
    #                              None,
    #                              tree)
    AddModelInstanceFromUrdfFile("resources/cartpole_disturbance.urdf", 
                                 FloatingBaseType.kFixed,
                                 None,
                                 tree)
    AddModelInstanceFromUrdfFile("resources/wall.urdf", 
                                 FloatingBaseType.kFixed,
                                 None,
                                 tree)

    plant = RigidBodyPlant(tree)
    context = plant.CreateDefaultContext()

    #d = 5.
    #d = 0.
    # d = 100
    # context.FixInputPort(0, BasicVector([d, 0.]))

    input_weight = 1.
    goal_weight = 30.
    #goal_weight = 5.
    #goal_x = 2.
    goal_x = 1.
    def running_cost(x, u, t):
        return input_weight * (u[0]**2 + u[1]**2) + goal_weight * (x[0] - goal_x)**2

    def final_cost(x, u, t):
        return 0.

    signed_dist_funcs = init_signed_dist_funcs()

    def final_state_constraint(x):
        pole_length = 0.5
        wall_x = 2.
        return x[0] + pole_length * sym.sin(x[1]) - wall_x

    #final_state_constraint = None

    planner = Planner(plant, context, running_cost, final_cost, signed_dist_funcs, final_state_constraint)

    #initial_state = (0., math.pi/2, 0., 0.)
    initial_state = (0., 0., 0., 0.)
    #final_state = (1., 3*math.pi/2., 0., 0.)
    final_state = None
    #final_state = (1.7, 0.64350111, 0.0, 0.0)
    #final_state = (1.5, 1.57079633, 0.0, 0.0)
    #duration_bounds = (3., 3.)
    #duration_bounds = (8., 8.)
    duration_bounds = (5., 5.)
    d = 0
    x_traj, total_cost = planner._solve_traj_opt(initial_state, final_state, duration_bounds, d, verbose=True)

    # from pydrake.all import PiecewisePolynomial
    # tmp = PiecewisePolynomial.FirstOrderHold([0., 1.],
    #                                          np.column_stack((initial_state,
    #                                                           initial_state)))
    # x_traj = x_traj + tmp

    duration = x_traj.end_time() - x_traj.start_time()
    print("Trajectory duration is {} s".format(duration))

    vis = PlanarRigidBodyVisualizer(tree, xlim=[-2.5, 2.5], ylim=[-1, 2.5])
    ani = vis.animate(x_traj, repeat=True)

    ani.save('out.mp4')

    plt.show()


def test_enumeration_planner():
    tree = RigidBodyTree()
    AddModelInstanceFromUrdfFile("resources/cartpole_disturbance.urdf", 
                                 FloatingBaseType.kFixed,
                                 None,
                                 tree)
    AddModelInstanceFromUrdfFile("resources/wall.urdf", 
                                 FloatingBaseType.kFixed,
                                 None,
                                 tree)

    plant = RigidBodyPlant(tree)
    context = plant.CreateDefaultContext()

    #d = 1.
    d = 2.
    #d = 5.
    #d = 10.
    #d = 0.

    '''
    input_weight = 1.
    goal_weight = 30.
    goal_x = 1.
    def running_cost(x, u, t):
        return input_weight * (u[0]**2 + u[1]**2) + goal_weight * (x[0] - goal_x)**2

    def final_cost(x, u, t):
        return 0.
    '''

    input_weight = 1.
    pole_upright_weight = 5.
    #goal_weight = 50.
    goal_weight = 30.
    goal_x = 1.

    def running_cost(x, u, t):
        return input_weight * (u[0]**2 + u[1]**2) + pole_upright_weight * x[1]**2

    def final_cost(x, u, t):
        return goal_weight * (x[0] - goal_x)**2

    def final_state_constraint(x):
        pole_length = 0.5
        wall_x = 2.
        return x[0] + pole_length * sym.sin(x[1]) - wall_x

    signed_dist_funcs = init_signed_dist_funcs()

    print("Disturbance is {} N".format(d))

    planner = Planner(plant, context, running_cost, final_cost, signed_dist_funcs, final_state_constraint)

    initial_state = (0., 0., 0., 0.)
    #tmin = 0.5
    tmin = 4.
    T = 8.
    #dt = 0.1
    #dc = 0.1
    dt = 0.5
    dc = 0.35
    x_traj_nc, x_traj_c = planner.plan(initial_state, tmin, T, dt, dc, d)

    vis = PlanarRigidBodyVisualizer(tree, xlim=[-2.5, 2.5], ylim=[-1, 2.5])
    ani = vis.animate(x_traj_nc, repeat=True)

    ani.save('disturbance_{}N.mp4'.format(d))

    plt.show()    


def main():
    #test_planner()
    test_enumeration_planner()


if __name__ == '__main__':
    main()
