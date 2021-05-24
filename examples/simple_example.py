# Import the package
from Buche_Silberstein_model_2020 import *

# Create the constitutive model
model = Buche_Silberstein_model_2020(number_of_links = 25, nondimensional_link_stiffness = 50)

# Define the applied deformation and specify the boundary conditions
def F(t): return 1 + t
F.type = 'uniaxial', '11'

# Solve the constitutive model (for the other deformations components, the stress, etc.)
model.solve(F, [0, 3])