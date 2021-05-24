# Buche_Silberstein_model_2020

[![GitHub release](https://mbuche.github.io/web/badges/releasev1.0.1.svg)](https://github.com/mbuche/Buche_Silberstein_model_2020/releases/) &nbsp; [![PyPI pyversions](https://mbuche.github.io/web/badges/python3.7.svg)](https://pypi.org/project/Buche-Silberstein-model-2020/) &nbsp; [![GitHub license](https://mbuche.github.io/web/badges/licenseMIT.svg)](https://github.com/mbuche/chain_breaking_polymer_networks/blob/master/LICENSE)

This is the Python package corresponding to: Buche, Michael R., and Meredith N. Silberstein. Statistical mechanical constitutive theory of polymer networks: The inextricable links between distribution, behavior, and ensemble. Physical Review E, 102, 012501 (2020).

[![DOI:10.1103/PhysRevE.102.012501](https://mbuche.github.io/web/badges/PhysRevE.102.012501.svg)](https://doi.org/10.1103/PhysRevE.102.012501) &nbsp; [![PhysRevE.102.012501](https://mbuche.github.io/web/badges/badgePRE102012501.svg)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.012501) &nbsp; [![arXiv:2004.07874](https://mbuche.github.io/web/badges/badgearXiv200407874.svg)](https://arxiv.org/abs/2004.07874)

# Installation

The package is on [PyPI](https://pypi.org/project/Buche-Silberstein-model-2020/) and can be installed using `pip`:

	pip install Buche_Silberstein_model_2020

It was written for `Python 3`, and uses `numpy` and `scipy`.

# Basic usage

The package is best imported using:

	from Buche_Silberstein_model_2020 import *
	
## Creating the constitutive model
	
The constitutive model for the network is created using:

	model = Buche_Silberstein_model_2020(number_of_links = 25, nondimensional_link_stiffness = 50)

The extensible freely-jointed chain (EFJC) model is automatically utilized. 
Note that the number of links in the chain (N_b) and the nondimensional link stiffness (kappa) have been specified as 25 and 50, respectively.
Also note the following optional keyword arguments:

* `shear_modulus`, which defaults to 1 (nondimensional stress);
* `bulk_modulus`, which defaults to 1000 times `shear_modulus` (effectively incompressible);
* `method`, which defaults to `'Gibbs_Legendre'` (the other option is `'Gaussian-Gibbs-Legendre'`, and `'Helmholtz'` will be implemented in a future release);
* `num_grid_romb`, which defaults to `1 + 2**9` and controls the number of grid points used for completing spatial integrals.

The constitutive model also inherits the following single-chain functions from the EFJC model, each as a function of the nondimensional extension (i.e. `model.eta(2.3)` gives the nondimensional force as the nondimensional extension 2.3):

* `eta`, the nondimensional force (Gibbs-Legendre method);
* `vartheta`, the nondimensional Helmholtz free energy per link (Gibbs-Legendre method);
* `P_eq`, the nondimensional equilibrium distribution (Gibbs-Legendre method);
* `g_eq`, the nondimensional equilibrium radial distribution function (Gibbs-Legendre method);
* `P_eq_Gaussian`, the nondimensional equilibrium distribution (Gibbs-Legendre method);
* `g_eq_Gaussian`, the nondimensional equilibrium radial distribution function (Gaussian-Gibbs-Legendre method).

## Deformation

After creating the constitutive model, we define the applied deformation and the associated traction boundary conditions:

	def F(t): return 1 + t
	F.type = 'uniaxial', '11'

This corresponds to uniaxial tension in the 1-direction. 
The solver will apply the boundary conditions and find the unknown deformation components. 
The solver also reduces spatial integration time (2D instead of 3D) through utilization of the symmetry conserved by the deformation in this case.

## Solution

The model is solved over a period of 3 seconds via

	model.solve(F, [0, 3])

An array of time values for the solution to be evaluated at can be given instead of a timespan.
The solution, attributed as `model.solution`, itself has several attributes (i.e. `model.solution.t` returns the attribute `t`):

* `t`, the times where the solution have been evaluated (shape N);
* `F`, the deformation gradient (shape 3 by 3 by N);
* `Cauchy_stress`, the true stress tensor (shape 3 by 3 by N);
* `nominal_stress`, the engineering stress tensor (shape 3 by 3 by N);
* `Hencky_strain`, the true strain tensor (shape 3 by 3 by N);

See [simple_example.py](examples/simple_example.py) for a self-contained version of this example.

# Parameter study

See [example_with_plots.py](examples/example_with_plots.py), where we generate (using `matplotlib`) Figures 3(a) and 4(a) of our [paper](https://doi.org/10.1103/PhysRevE.102.012501).

<table>
	<tr>
		<td> 
			<img src="https://github.com/mbuche/Buche_Silberstein_model_2020/blob/main/examples/g_eq.png" alt="g_eq" style="width: 250px;"/>
		</td>
		<td> 
			<img src="https://github.com/mbuche/Buche_Silberstein_model_2020/blob/main/examples/stress.png" alt="stress" style="width: 250px;"/>
		</td>
	</tr>
</table>
