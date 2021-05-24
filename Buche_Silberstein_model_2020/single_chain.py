################################################################################################################################
# General setup
################################################################################################################################

# Import miscellaneous functions and variables
from Buche_Silberstein_model_2020.miscellaneous import *

# Import libraries
import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar

################################################################################################################################
# General single-chain class
################################################################################################################################

class single_chain:

	# Initialization
	def __init__(self, **kwargs):

		# Default parameter values
		self.epsabs = 1e-3
		self.epsrel = 1e-3
		self.num_interp = int(3e3)
		self.interp_kind_1D = 'cubic'
		self.cutoff_for_log_over_sinh = 3e1
		self.cutoff_stretch_for_harmonic = 3
		self.N_b = 88
		self.kappa = 888

		# Retrieve specified parameters
		for key, value in kwargs.items():
			exec("self.%s = value" % key, {'self': self, 'value': value})

	# Function to invert a function
	def inv_fun_1D(self, x_query, fun, bounds = None):
		return root_scalar(lambda x: fun(x) - x_query, x0 = x_query*0.95, x1 = x_query/0.95).root

	# Function to avoid overflow when computing ln(x/sinh(x))
	def log_over_sinh(self, x):

		# Determine when argument is sufficiently large
		where_x_large = np.nan_to_num(x, nan = -1) > self.cutoff_for_log_over_sinh
		log_of_x_over_sinh_x = np.zeros(x.shape)

		# Use asymptotic relation valid for sufficiently large arguments
		if where_x_large.any():
			log_of_x_over_sinh_x[where_x_large] = np.log(2*x[where_x_large]) - x[where_x_large]

		# Compute analytically otherwise, and zero where argument is zero
		where_x_zero = x == 0
		where_compute = ~(where_x_large + where_x_zero)
		if where_compute.any():
			log_of_x_over_sinh_x[where_compute] = np.log(x[where_compute]/np.sinh(x[where_compute]))
		return log_of_x_over_sinh_x

	# Hyperbolic cotangent function
	def coth_safe(self, eta):
		eta = np.where(eta == 0, minimum_float, eta)
		return 1/np.tanh(eta)

	# Langevin function
	def Langevin(self, eta):
		eta = np.where(eta == 0, minimum_float, eta)
		return 1/np.tanh(eta) - 1/eta

	# Inverse Langevin function
	def inv_Langevin(self, gamma):
		return self.inv_fun_1D(gamma, self.Langevin)

################################################################################################################################
# Extensible freely-joined chain model
################################################################################################################################

class EFJC(single_chain):

	# For more information, see:
	# 	Analytical results of the extensible freely jointed chain model
	# 	Alessandro Fiasconaro and Fernando Falo
	#	Physica A, 532, 121929 (2019)
	# 	doi.org/10.1016/j.physa.2019.121929
	# See also:
	#	Statistical mechanical constitutive theory of polymer networks: 
	#		The inextricable links between distribution, behavior, and ensemble
	# 	Michael R. Buche and Meredith N. Silberstein
	# 	Physical Review E, 102, 012501 (2020)
	# 	doi.org/10.1103/PhysRevE.102.012501

	# Initialization
	def __init__(self, **kwargs):

		# Retrieve default and specified parameter values
		single_chain.__init__(self, **kwargs)

		# Nondimensional mechanical response of a single chain
		def gamma_fun(eta):
			coth = self.coth_safe(eta)
			L = self.Langevin(eta)
			return L + eta/self.kappa*(1 + (1 - L*coth)/(1 + eta/self.kappa*coth))

		# Compute and store the inverted nondimensional mechanical response to interpolate from
		self.gamma_store = np.linspace(0, self.cutoff_stretch_for_harmonic, self.num_interp)
		self.eta_store = np.zeros(self.gamma_store.size)
		for i in range(1, len(self.gamma_store)):
			self.eta_store[i] = self.inv_fun_1D(self.gamma_store[i], gamma_fun)

		# Function to interpolate from the inverted nondimensional mechamical response of the chain
		def eta_fun(gamma):
			if isinstance(gamma, np.ndarray):
				eta_out = np.zeros(gamma.shape)
				harmonic_region = gamma > self.cutoff_stretch_for_harmonic
				eta_out[harmonic_region] = self.kappa*(gamma[harmonic_region] - 1)
				eta_out[~harmonic_region] = np.interp(gamma[~harmonic_region], self.gamma_store, self.eta_store)
			else:
				if gamma > self.cutoff_stretch_for_harmonic:
					eta_out = self.kappa*(gamma - 1)
				else:
					eta_out = np.interp(gamma, self.gamma_store, self.eta_store)
			return eta_out

		# Nondimensional Helmholtz free energy per link of a single chain
		def vartheta_fun(gamma):

			# Compute mechanical response
			eta = np.array(eta_fun(gamma))
			eta[eta == 0] = minimum_float

			# Compute nondimensional Helmholtz free energy per link
			coth = self.coth_safe(eta)
			L = self.Langevin(eta)
			return eta*L + self.log_over_sinh(eta) - np.log(1 + eta/self.kappa*coth) \
				+ eta**2/self.kappa/2*(1/2 + (1 - L*coth)/(1 + eta/self.kappa*coth))

		# Nondimensional equilibrium distribution function
		def P_eq_fun(gamma, normalization = 1):

			# Compute nondimensional Helmholtz free energy per link
			vartheta = vartheta_fun(gamma)

			# Compute and return P_eq
			return np.exp(-self.N_b*vartheta)/normalization

		# Nondimensional equilibrium radial distribution function
		def g_eq_fun(gamma, normalization = 1):
			return 4*np.pi*gamma**2*P_eq_fun(gamma, normalization)

		# Normalize the equilibrium distribution
		P_eq_normalization = quad(g_eq_fun, 0, np.inf, epsabs = self.epsabs, epsrel = self.epsrel, full_output = 1)[0]

		# Compute average nondimensional end-to-end length at equilibrium
		integrand = lambda gamma: gamma*g_eq_fun(gamma, normalization = P_eq_normalization)
		self.average_gamma_eq = quad(integrand, 0, np.inf, epsabs = self.epsabs, epsrel = self.epsrel, full_output = 1)[0]

		# Nondimensional equilibrium (Gaussian) distribution function valid in the limit N_b -> inf
		def P_eq_Gaussian_fun(gamma):
			c_kappa = self.kappa*(self.kappa + 1)/(self.kappa*self.kappa + 6*self.kappa + 3)
			return (3*c_kappa*self.N_b/(2*np.pi))**(3/2)*np.exp(-3/2*c_kappa*self.N_b*gamma*gamma)

		# Nondimensional equilibrium radial (Gaussian) distribution function valid in the limit N_b -> inf
		def g_eq_Gaussian_fun(gamma, normalization = 1):
			return 4*np.pi*gamma**2*P_eq_Gaussian_fun(gamma)

		# Retain each single-chain function
		self.eta = eta_fun
		self.vartheta = vartheta_fun
		self.P_eq = lambda gamma: P_eq_fun(gamma, normalization = P_eq_normalization)
		self.g_eq = lambda gamma: g_eq_fun(gamma, normalization = P_eq_normalization)
		self.P_eq_Gaussian = lambda gamma: P_eq_Gaussian_fun(gamma)
		self.g_eq_Gaussian = lambda gamma: g_eq_Gaussian_fun(gamma)
