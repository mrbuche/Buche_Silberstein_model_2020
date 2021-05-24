################################################################################################################################
# Miscellaneous functions and variables
################################################################################################################################

# Import libraries
import sys
import numpy as np
import scipy.linalg as la

# Use to avoid overflow
minimum_float = sys.float_info.min
maximum_float = sys.float_info.max
minimum_exponent = np.log(minimum_float)/np.log(10)
maximum_exponent = np.log(maximum_float)/np.log(10)

# The identity tensors
identity_tensor = np.diag([1, 1, 1])
fourth_rank_identity_tensor = np.tensordot(identity_tensor, identity_tensor, axes = 0)

# Function returning the deviatoric portion of a tensor
def dev(A):
	return A - A.trace()/3*identity_tensor

# Function to check symmetries of deformation gradient
def symmetry_check(F):
	diagonal_check = np.count_nonzero(F - np.diag(np.diagonal(F)))
	if diagonal_check == 0:
		if np.isclose(F[1, 1]**2, 1/F[0, 0]) and np.isclose(F[1, 1], F[2, 2]):
			return 'uniaxial'
		elif np.isclose(F[2, 2], 1/F[0, 0]**2) and np.isclose(F[1, 1], F[0, 0]):
			return 'equibiaxial'
		else:
			return 'diagonal'
	else:
		return 'none'
