################################################################################################################################
# General setup
################################################################################################################################

# Import single-chain model library
from Buche_Silberstein_model_2020.single_chain import *

# Import libraries
import numpy as np
import scipy.linalg as la
from scipy.integrate import romb
from scipy.optimize import root_scalar

################################################################################################################################
# Solve a model with no internal state variables
################################################################################################################################

class solver:

	# Currently only implemented in F principal coordinates 
	# Also assumes initially isotropic and some symmetry is preserved by the deformation
	def solve(self, F_applied, t_span, *args, **kwargs):

		# Residual function to apply traction boundary conditions
		if F_applied.type[0] == 'uniaxial':
			def guess(t):
				return 1/np.sqrt(F_applied(t))
			app_comp = int(F_applied.type[1][0]) - 1
			def residual(t, F_traction_free, return_F = False):
				F = np.zeros((3, 3))
				F[app_comp, app_comp] = F_applied(t)
				F[app_comp - 1, app_comp - 1] = F_traction_free
				F[app_comp - 2, app_comp - 2] = F_traction_free
				if return_F is True:
					return F
				else:
					return self.Cauchy_stress(t, F)[app_comp - 1, app_comp - 1]

		# Default to a certain number of points in time
		if len(t_span) == 2:
			t_span = np.linspace(t_span[0], t_span[1], 23)
		class solution:
			pass
		solution.t = t_span

		# Function to solve for F components using boundary conditions
		def F_fun(t):
			F_traction_free = root_scalar(lambda F: residual(t, F), x0 = guess(t)*0.95, x1 = guess(t)/0.95).root
			return residual(t, F_traction_free, return_F = True)

		# Recompute and store all components of deformation gradient 
		solution.F = np.zeros((3, 3, len(solution.t)))
		for index_t in range(len(solution.t)):
			solution.F[:, :, index_t] = F_fun(solution.t[index_t])

		# Compute relevant stresses
		solution.Cauchy_stress = np.zeros((3, 3, len(solution.t)))
		solution.nominal_stress = np.zeros((3, 3, len(solution.t)))
		solution.Hencky_strain = np.zeros((3, 3, len(solution.t)))
		for index_t in range(len(solution.t)):
			solution.Cauchy_stress[:, :, index_t] = \
				self.Cauchy_stress(solution.t[index_t], solution.F[:, :, index_t])
			solution.nominal_stress[:, :, index_t] = \
				la.det(solution.F[:, :, index_t])*la.inv(solution.F[:, :, index_t]).dot(solution.Cauchy_stress[:, :, index_t])
			solution.Hencky_strain[:, :, index_t] = \
				la.logm(la.sqrtm(solution.F[:, :, index_t].dot(solution.F[:, :, index_t].T)))

		# Return solution object with t, F, and stresses as attributes
		self.solution = solution

################################################################################################################################
# General class for elastic models
################################################################################################################################

class elastic(solver):

	# Initialization
	def __init__(self, **kwargs):

		# Default parameter values
		self.shear_modulus = 1
		self.number_of_links = 88
		self.nondimensional_link_stiffness = 888
		self.num_grid_romb = 1 + 2**9
		self.method = 'Gibbs-Legendre'

		# Retrieve specified parameters
		for key, value in kwargs.items():
			exec("self.%s = value" % key, {'self': self, 'value': value})

	# Check if bulk modulus specified, otherwise ensure much larger than modulus
	def bulk_modulus_fun(self):
		try:
			self.bulk_modulus
		except AttributeError:
			return 1e3*self.shear_modulus
		else:
			return self.bulk_modulus

	# Default spherical term in Cauchy stress for effectively incompressible models
	def spherical_Cauchy_stress(self, F):
		J = la.det(F)
		return self.bulk_modulus_fun()*(J - 1)*identity_tensor

################################################################################################################################
# Neo-Hookean model
################################################################################################################################

class Neo_Hookean(elastic):

	# For more information, see:
	#	Large elastic deformations of isotropic materials. I. Fundamental concepts.
	# 	Ronald S. Rivlin
	# 	Philosophical Transactions of the Royal Society of London. 
	#		Series A, Mathematical and Physical Sciences, 240, 822 (1948)
	# 	doi.org/10.1098/rsta.1948.0002

	# Initialization
	def __init__(self, **kwargs):

		# Retrieve default and specified parameter values
		elastic.__init__(self, **kwargs)

	# Stress due to the applied deformation
	def Cauchy_stress(self, *args):

		# Kinematics
		F = args[-1]
		J = la.det(F)
		B = F.dot(F.T)
		dev_B_bar = dev(B)/J**(2/3)

		# Constitutive relation for the stress
		return self.shear_modulus*dev_B_bar/J + self.spherical_Cauchy_stress(F)

################################################################################################################################
# Buche-Silberstein model
################################################################################################################################

class Buche_Silberstein_model_2020(elastic, EFJC):

	# For more information, see:
	#	Statistical mechanical constitutive theory of polymer networks: 
	#		The inextricable links between distribution, behavior, and ensemble
	# 	Michael R. Buche and Meredith N. Silberstein
	# 	Physical Review E, 102, 012501 (2020)
	# 	doi.org/10.1103/PhysRevE.102.012501

	# Initialization
	def __init__(self, **kwargs):

		# Retrieve default and specified parameter values
		elastic.__init__(self, **kwargs)
		EFJC.__init__(self, N_b = self.number_of_links, kappa = self.nondimensional_link_stiffness)

		# Mark grid not initialized
		self.grid_initialized = False

		# Will implement the Helmholtz method here in the future
		if self.method == 'Helmholtz':
			pass
		elif self.method == 'Gibbs-Legendre':
			self.P_eq_used = self.P_eq
		elif self.method == 'Gaussian-Gibbs-Legendre':
			self.P_eq_used = self.P_eq_Gaussian

	def initialized_grid(self, symmetry):

		# Mark grid as initialized using the given symmetry
		self.grid_initialized = symmetry
		
		if symmetry == 'uniaxial':
			w = np.linspace(0, 1, 1 + self.num_grid_romb)[:-1]
			self.dw = w[1] - w[0]
			W_z, W_r = np.meshgrid(w, w)
			self.Z = np.arctanh(W_z)
			self.R = np.arctanh(W_r)
			GAMMA = np.sqrt(self.Z*self.Z + self.R*self.R)
			ELEMENT = 4*np.pi*self.R/(1 - W_z*W_z)/(1 - W_r*W_r)
			GEOM = np.zeros((self.num_grid_romb, self.num_grid_romb, 3, 3))
			GEOM[:, :, 0, 0] = self.Z*self.Z
			GEOM[:, :, 1, 1] = self.R*self.R/2
			GEOM[:, :, 2, 2] = self.R*self.R/2
			ETA_OVER_GAMMA = np.zeros(GAMMA.shape)
			ETA_OVER_GAMMA[GAMMA != 0] = self.eta(GAMMA[GAMMA != 0])/GAMMA[GAMMA != 0]
			self.ELEMENT_STRESS = (self.shear_modulus*self.number_of_links*ETA_OVER_GAMMA*ELEMENT)[:, :, None, None]*GEOM

		# No (None) symmetry
		else:
			w = np.linspace(-1, 1, 2 + self.num_grid_romb)[1:-1]
			self.dw = w[1] - w[0]
			W_x, W_y, W_z = np.meshgrid(w, w, w)
			self.X = np.arctanh(W_x)
			self.Y = np.arctanh(W_y)
			self.Z = np.arctanh(W_z)
			GAMMA = np.sqrt(self.X*self.X + self.Y*self.Y + self.Z*self.Z)
			ELEMENT = 1/(1 - W_x*W_x)/(1 - W_y*W_y)/(1 - W_z*W_z)
			GEOM = np.zeros((self.num_grid_romb, self.num_grid_romb, self.num_grid_romb, 3, 3))
			GEOM[:, :, :, 0, 0] = self.X*self.X
			GEOM[:, :, :, 0, 1] = self.X*self.Y
			GEOM[:, :, :, 0, 2] = self.X*self.Z
			GEOM[:, :, :, 1, 0] = self.Y*self.X
			GEOM[:, :, :, 1, 1] = self.Y*self.Y
			GEOM[:, :, :, 1, 2] = self.Y*self.Z
			GEOM[:, :, :, 2, 0] = self.Z*self.X
			GEOM[:, :, :, 2, 1] = self.Z*self.Y
			GEOM[:, :, :, 2, 2] = self.Z*self.Z
			ETA_OVER_GAMMA = np.zeros(GAMMA.shape)
			ETA_OVER_GAMMA[GAMMA != 0] = self.eta(GAMMA[GAMMA != 0])/GAMMA[GAMMA != 0]
			self.ELEMENT_STRESS = (self.shear_modulus*self.number_of_links*ETA_OVER_GAMMA*ELEMENT)[:, :, :, None, None]*GEOM

	# Stress due to the applied deformation
	def Cauchy_stress(self, *args):

		# Kinematics
		F = args[-1]
		J = la.det(F)
		F_bar = F/J**(1/3)

		if symmetry_check(F_bar) == 'uniaxial':

			# Only initialize the grid if it has not been initialized already using this symmetry
			if self.grid_initialized != 'uniaxial':
				self.initialized_grid('uniaxial')

			# Relatively-deformed coordinates
			GAMMA_F = np.sqrt((self.Z/F_bar[0, 0])**2 + F_bar[0, 0]*self.R*self.R)

			# Integrate for and return the stress
			INTEGRAND = self.P_eq_used(GAMMA_F)[:, :, None, None]*self.ELEMENT_STRESS
			return self.spherical_Cauchy_stress(F) + romb(romb(INTEGRAND, dx = self.dw, axis = 0), dx = self.dw, axis = 0)/J

		# General, no symmetry utilized
		else:

			# Only initialize the grid if it has not been initialized already using this symmetry
			if self.grid_initialized != None:
				self.initialized_grid(None)

			# Relatively-deformed coordinates
			F_bar_inv = la.inv(F_bar)
			X_F = F_bar_inv[0, 0]*self.X + F_bar_inv[0, 1]*self.Y + F_bar_inv[0, 2]*self.Z
			Y_F = F_bar_inv[1, 0]*self.X + F_bar_inv[1, 1]*self.Y + F_bar_inv[1, 2]*self.Z
			Z_F = F_bar_inv[2, 0]*self.X + F_bar_inv[2, 1]*self.Y + F_bar_inv[2, 2]*self.Z
			GAMMA_F = np.sqrt(X_F*X_F + Y_F*Y_F + Z_F*Z_F)

			# Integrate for and return the stress
			INTEGRAND = self.P_eq_used(GAMMA_F)[:, :, :, None, None]*self.ELEMENT_STRESS
			return self.spherical_Cauchy_stress(F) + \
				romb(romb(romb(INTEGRAND, dx = self.dw, axis = 0), dx = self.dw, axis = 0), dx = self.dw, axis = 0)/J
