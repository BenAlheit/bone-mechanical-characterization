import numpy as np
import matplotlib.pyplot as plt
import config.plotting_config
from viscoelastic_framework_1D import LinearVEMaterial, LinearEMOdWVMaterial, BackwardEuler, VE1D, Lagrange
from linear_analytical_constant_rate import sigma as sigma_an

E_inf = 2e9
E_s = 1e9
tau = 1e-3
n = 0.1
K = E_s * tau
eps_dot = np.logspace(-4, 6, 1000)
eps = np.linspace(0.001, 0.01, 5)

# epss = np.linspace(0, 0.01, 500)
# epss = np.linspace(0, 0.01, 50)
epss = np.linspace(0, 0.01, 100)
eps_dots = np.logspace(-3, 5, 50)
E_ss = np.linspace(0.8e9, 1.2e9, 5)
taus = np.logspace(-5, -1, 5)

ns = np.logspace(-1, 1, 5)

fig_path = '../../written-work/viscoelasticity/figures/'


def sigma(E_inf, E_s, K, n, eps, eps_dot):
    ts = eps / eps_dot
    sigs = E_inf * eps
    eds = np.zeros(eps.size)

    mat = LinearEMOdWVMaterial(E_s, K, n)
    # ta = BackwardEuler()
    ta = Lagrange(2)
    ve1d = VE1D(mat, ta)

    for it in range(1, ts.size):
        ve1d.increment_model(ts[it], eps[it])
        sigs[it] += ve1d.sig_e
        eds[it] = ve1d.ta.a()

    return sigs, eds


def format_and_save(path, x_label=r'$\dot{\varepsilon}$', save=False):
    plt.grid()
    plt.legend()
    plt.xlabel(x_label)
    plt.ylabel(r'$\sigma$')
    plt.tight_layout()
    if save:
        plt.savefig(path)


for i_n in ns:
    sigmas = []
    for eps_dot in eps_dots:
        print(eps_dot)
        sigmas.append(sigma(E_inf, E_s, K, i_n, epss, eps_dot)[0][-1])

    plt.semilogx(eps_dots, sigmas, label=r'$n = ' + str(i_n) + '$')

format_and_save(f'{fig_path}/modw.pdf', save=True)

# plt.plot(epss, sigma_an(E_inf, E_s, tau, epss, eps_dots[-1]), label='Analytical', )
# plt.plot(epss, sigma(E_inf, E_s, tau, epss, eps_dots[-1]), label='Framework')
#
# format_and_save(f'{fig_path}e-stress.pdf', r'$\varepsilon$')

# for ep in eps:
#     plt.semilogx(eps_dot, sigma(E_inf, E_s, tau, ep, eps_dot), label=r'$\varepsilon = ' + str(ep) + '$')
#
# format_and_save(f'{fig_path}e-dot-stress.pdf')
# plt.figure()
#
# for i_tau in taus:
#     plt.semilogx(eps_dot, sigma(E_inf, E_s, i_tau, eps[-1], eps_dot), label=r'$\tau = ' + str(i_tau) + '\,s$')
#
# format_and_save(f'{fig_path}e-dot-stress-tau.pdf')
# plt.figure()
#
# for i_E_s in E_ss:
#     plt.semilogx(eps_dot, sigma(E_inf, i_E_s, tau, eps[-1], eps_dot), label=r'$E_s = ' + str(i_E_s) + '\,s$')
#
# format_and_save(f'{fig_path}e-dot-stress-E.pdf')
# plt.figure()

# for i_eps_dot in eps_dots:
#     sigs, eds = sigma(E_inf, E_s, K, n, epss, i_eps_dot)
#     # plt.plot(epss, eds, label=r'$\dot{\varepsilon} = ' + str(i_eps_dot) + '\,/s$')
#     plt.plot(epss, sigma_an(E_inf, E_s, tau, epss, i_eps_dot), label='Analytical', color='blue')
#     plt.plot(epss, sigs, label='Framework', color='red', linestyle='--')
#
# format_and_save(f'{fig_path}e-stress.pdf', r'$\varepsilon$')

plt.show()
