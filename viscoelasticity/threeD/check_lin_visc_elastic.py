import solid_mechanics as sm
import numpy as np
import matplotlib.pyplot as plt
import config.plotting_config as plt_cfig

E_inf = 2e9
v_inf = 0.001
E_s = 1e9
v_s = 0.001
tau = 1e-3

eps_dot = np.logspace(-4, 6, 1000)
eps = np.linspace(0.001, 0.01, 5)

epss = np.linspace(0, 0.01, 200)
eps_dots = np.logspace(-3, 5, 1000)


fig_path = '../../written-work/viscoelasticity/scratch/figures/'


def e_v_2_lam_mu(E, v):
    return E*v/((1+v)*(1-2*v)), E/(2*(1+v))

lambda_inf, mu_inf = e_v_2_lam_mu(E_inf, v_inf)
lam_s, mu_s = e_v_2_lam_mu(E_s, v_s)

eta_d = E_s * tau


def sigma(lambda_inf, mu_inf, eta_d, mu_s, eps, eps_dot):
    ts = eps / eps_dot
    sigs = np.zeros(eps.shape)
    eds = np.zeros(eps.size)

    # mat = LinearECarreuVMaterial(E_s, eta_0, eta_inf, lam, n, a)
    # ta = BackwardEuler()
    # ta = Lagrange(2)
    # ve1d = VE1D(mat, ta)

    eps_n1 = np.zeros(6)
    eps_dn = np.zeros(6)

    for it in range(1, ts.size):
        dt = ts[it] - ts[it-1]
        eps_n1[0] = eps[it]
        eps_dn1, sigma = sm.solid_mechanics.increment_lin_visc_elastic(mu_inf, lambda_inf, eta_d, mu_s, eps_n1, eps_dn, dt)
        sigs[it] = sigma[0]
        eds[it] = eps_dn1[0]

        eps_dn = eps_dn1.copy()

    return sigs, eds



sigmas = []
for eps_dot in eps_dots:
    print(eps_dot)
    sigmas.append(sigma(lambda_inf, mu_inf, eta_d, mu_s, epss, eps_dot)[0][-1])

# plt.semilogx(eps_dots, sigmas, label=r'$n = ' + str(i_n) + '$')
plt.semilogx(eps_dots, sigmas)

plt_cfig.format_and_save(f'{fig_path}/odw.pdf', save=False)

plt.show()