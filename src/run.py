import math
import numpy as np
import matplotlib.pyplot as plt

from pydrake.examples.pendulum import (PendulumPlant, PendulumState)
from pydrake.all import (DirectCollocation, PiecewisePolynomial,
                         SolutionResult)
from visualizer import PendulumVisualizer

if __name__ == "__main__":
	# planner = Planner()
