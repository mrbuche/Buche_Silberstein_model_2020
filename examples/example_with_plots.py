# Import the package
from Buche_Silberstein_model_2020 import *

# Import matplotlib for plotting and saving
import matplotlib.pyplot as plt
def save_current_figure(xlabel, ylabel, name):
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.tight_layout()
	plt.show()
	fig = plt.gcf()
	plt.savefig(name)
	plt.close()

# Compare equilibrium radial distribution functions as the number of links increases
number_of_links_list = [3, 5, 10, 25]
gamma_plotting = np.linspace(0, 1.6, 1000)
for number_of_links in number_of_links_list:
	model = Buche_Silberstein_model_2020(number_of_links = number_of_links, nondimensional_link_stiffness = 50)
	plt.plot(gamma_plotting, model.g_eq(gamma_plotting), color = 'blue')
	plt.plot(gamma_plotting, model.g_eq_Gaussian(gamma_plotting), color = 'red')
plt.legend(['Gibbs-Legendre', 'Gaussian'])
save_current_figure('Nondimensional single-chain extension', 'Nondimensional radial distribution function', 'g_eq.png')

# Compare uniaxial tension results as the number of links increases
def F(t): return 1 + t
F.type = 'uniaxial', '11'
model_0 = Neo_Hookean()
model_0.solve(F, [0, 3])
for number_of_links in number_of_links_list[1:]:
	model_1 = Buche_Silberstein_model_2020(number_of_links = number_of_links, nondimensional_link_stiffness = 50)
	model_2 = Buche_Silberstein_model_2020(number_of_links = number_of_links, nondimensional_link_stiffness = 50, \
		method = 'Gaussian-Gibbs-Legendre')
	model_1.solve(F, [0, 3])
	model_2.solve(F, [0, 3])
	plt.plot(model_1.solution.F[0, 0, :], model_1.solution.Cauchy_stress[0, 0, :], color = 'blue')
	plt.plot(model_2.solution.F[0, 0, :], model_2.solution.Cauchy_stress[0, 0, :], color = 'red')
	if number_of_links == number_of_links_list[1]:
		plt.plot(model_0.solution.F[0, 0, :], model_0.solution.Cauchy_stress[0, 0, :], color = 'green')
plt.xlim([1, 4])
plt.ylim([0, 35])
plt.legend(['Gibbs-Legendre', 'Gaussian-Gibbs-Legendre', 'Neo-Hookean'])
save_current_figure('Applied stretch (uniaxial tension)', 'Nondimensional stress', 'stress.png')
