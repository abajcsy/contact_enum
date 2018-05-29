import math
import numpy as np

class Planner:

	def enumerate_modes(self, contact_surf, dt, dc, T):
		return None

	def plan_hybrid_traj(self, xinit, d, J, c, t, dt):
		return None

	def plan(self)
		# C = self.enumerate_modes()
		# for (c,t) in C
		# 	S = S union {self.plan_hybrid_traj()}
		# return argmin(J, (x,u))