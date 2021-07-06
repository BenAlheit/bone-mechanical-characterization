import numpy as np
import matplotlib.pyplot as plt
import config.plotting_config
from viscoelasticity.oneD.viscoelastic_framework_1D import LinearECarreuVMaterial, BackwardEuler, VE1D, Lagrange
from viscoelasticity.oneD.linear_analytical_constant_rate import sigma as sigma_an

E_inf = 2e9
E_s = 1e9
tau = 1e-7
n = 0.33
eta_0 = E_s * tau
eta_inf = eta_0*1e5
a = 1
# a = 2
# lam = 0.1
lam = 1
save = True

eps_dot = np.logspace(-4, 6, 1000)
eps = np.linspace(0.001, 0.01, 5)

# epss = np.linspace(0, 0.01, 500)
epss = np.linspace(0, 0.01, 100)
eps_dots = np.logspace(-3, 5, 50)
E_ss = np.linspace(0.8e9, 1.2e9, 5)
taus = np.logspace(-5, -1, 5)
a_s = np.logspace(-1, 2, 5)

# ns = np.logspace(-5, 0, 5)
ns = np.linspace(0.98, 0.9999, 5)
lams = np.logspace(-1, 1, 5)
eta_inf_mults = np.logspace(3, 10, 5)

fig_path = '../../written-work/viscoelasticity/scratch/figures/'


def sigma(E_inf, E_s, eta_0, eta_inf, lam, n, a, eps, eps_dot):
    ts = eps / eps_dot
    sigs = E_inf * eps
    eds = np.zeros(eps.size)

    mat = LinearECarreuVMaterial(E_s, eta_0, eta_inf, lam, n, a)
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


for i_lam in lams:
    sigmas = []
    for eps_dot in eps_dots:
        print(eps_dot)
        sigmas.append(sigma(E_inf, E_s, eta_0, eta_inf, i_lam, n, a, epss, eps_dot)[0][-1])

    plt.semilogx(eps_dots, sigmas, label=r'$\lambda = ' + str(i_lam) + '\,s$')

format_and_save(f'{fig_path}/carreau-lams.pdf', save=save)

plt.figure()

for i_eta_inf_mult in eta_inf_mults:
    sigmas = []
    for eps_dot in eps_dots:
        print(eps_dot)
        sigmas.append(sigma(E_inf, E_s, eta_0, i_eta_inf_mult * eta_0, lam, n, a, epss, eps_dot)[0][-1])

    plt.semilogx(eps_dots, sigmas, label=r'$\eta^\infty/\eta_0 = ' + str(i_eta_inf_mult) + '$')

format_and_save(f'{fig_path}/carreau-etas.pdf', save=save)
plt.figure()

for i_n in ns:
    sigmas = []
    for eps_dot in eps_dots:
        print(eps_dot)
        sigmas.append(sigma(E_inf, E_s, eta_0, eta_inf, lam, i_n, a, epss, eps_dot)[0][-1])

    plt.semilogx(eps_dots, sigmas, label=r'$n = ' + str(i_n) + '$')

format_and_save(f'{fig_path}/carreau-ns.pdf', save=save)

plt.figure()

for i_a in a_s:
    sigmas = []
    for eps_dot in eps_dots:
        print(eps_dot)
        sigmas.append(sigma(E_inf, E_s, eta_0, eta_inf, lam, n, i_a, epss, eps_dot)[0][-1])

    plt.semilogx(eps_dots, sigmas, label=r'$a = ' + str(i_a) + '$')

format_and_save(f'{fig_path}/carreau-as.pdf', save=save)


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
