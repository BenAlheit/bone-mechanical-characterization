import solid_mechanics as sm
import numpy as np
import matplotlib.pyplot as plt
import config.plotting_config as plt_cfig

E_inf = 2e9
v_inf = 0.001
E_s = 1e9
v_s = 0.001
tau = 1e-3

# # eps_dot = np.logspace(-6, 2, 1000)
# eps = np.linspace(0.001, 0.01, 5)
#
# epss = np.linspace(0, 0.01, 500)
# # eps_dots = np.logspace(-4, 2, 1000)

save=False

fig_path = '../../written-work/viscoelasticity/scratch/figures/'


def e_v_2_lam_mu(E, v):
    return E*v/((1+v)*(1-2*v)), E/(2*(1+v))

lambda_inf, mu_inf = e_v_2_lam_mu(E_inf, v_inf)
lam_s, mu_s = e_v_2_lam_mu(E_s, v_s)

eta_d = E_s * tau

n = 0.33
# n = 1
eta_0 = E_s * tau
# eta_inf = eta_0
eta_inf = eta_0*1e5
a = 1

lam = 1

epss = np.linspace(0, 0.01, 200)
eps_dots = np.logspace(-5, 2, 150)
E_ss = np.linspace(0.8e9, 1.2e9, 5)
taus = np.logspace(-5, -1, 5)
a_s = np.logspace(-1, 1.5, 5)

# ns = np.logspace(-5, 0, 5)
ns = np.linspace(0.98, 0.9999, 5)
lams = np.logspace(-1, 1, 5)
eta_inf_mults = np.logspace(3, 10, 5)


def sigma(lambda_inf, mu_inf, eta_0, eta_inf, a, lam, n, eps, eps_dot):
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
        # eps_dn1, sigma = sm.solid_mechanics.increment_lin_visc_elastic(mu_inf, lambda_inf, eta_d, mu_s, eps_n1, eps_dn, dt)
        eps_dn1, sigma = sm.solid_mechanics.increment_carraeu_visc_elastic(mu_inf,
                                                                           lambda_inf,
                                                                           lam,
                                                                           n, a,
                                                                           eta_0, eta_inf,
                                                                           mu_s, eps_n1, eps_dn, dt)
        sigs[it] = sigma[0]
        eds[it] = eps_dn1[0]

        eps_dn = eps_dn1.copy()

    return sigs, eds



sigmas = []
for eps_dot in eps_dots:
    print(eps_dot)
    sigmas.append(sigma(lambda_inf, mu_inf, eta_0, eta_inf, a, lam, n, epss, eps_dot)[0][-1])

# plt.semilogx(eps_dots, sigmas, label=r'$n = ' + str(i_n) + '$')
plt.semilogx(eps_dots, sigmas)

# plt_cfig.format_and_save(f'{fig_path}/odw.pdf', save=False)


for i_lam in lams:
    sigmas = []
    for eps_dot in eps_dots:
        print(eps_dot)
        # sigmas.append(sigma(E_inf, E_s, eta_0, eta_inf, i_lam, n, a, epss, eps_dot)[0][-1])
        sigmas.append(sigma(lambda_inf, mu_inf, eta_0, eta_inf, a, i_lam, n, epss, eps_dot)[0][-1])


    plt.semilogx(eps_dots, sigmas, label=r'$\lambda = ' + str(i_lam) + '\,s$')

plt_cfig.format_and_save(f'{fig_path}/carreau-3d-lams.pdf', save=save)

plt.figure()

for i_eta_inf_mult in eta_inf_mults:
    sigmas = []
    for eps_dot in eps_dots:
        print(eps_dot)
        # sigmas.append(sigma(E_inf, E_s, eta_0, i_eta_inf_mult * eta_0, lam, n, a, epss, eps_dot)[0][-1])
        sigmas.append(sigma(lambda_inf, mu_inf, eta_0,  i_eta_inf_mult * eta_0, a, lam, n, epss, eps_dot)[0][-1])


    plt.semilogx(eps_dots, sigmas, label=r'$\eta^\infty/\eta_0 = ' + str(i_eta_inf_mult) + '$')

plt_cfig.format_and_save(f'{fig_path}/carreau-3d-etas.pdf', save=save)
plt.figure()

for i_n in ns:
    sigmas = []
    for eps_dot in eps_dots:
        print(eps_dot)
        # sigmas.append(sigma(E_inf, E_s, eta_0, eta_inf, lam, i_n, a, epss, eps_dot)[0][-1])
        sigmas.append(sigma(lambda_inf, mu_inf, eta_0,  eta_inf, a, lam, i_n, epss, eps_dot)[0][-1])

    plt.semilogx(eps_dots, sigmas, label=r'$n = ' + str(i_n) + '$')



plt_cfig.format_and_save(f'{fig_path}/carreau-3d-ns.pdf', save=save)

plt.figure()

for i_a in a_s:
    sigmas = []
    for eps_dot in eps_dots:
        print(eps_dot)
        # sigmas.append(sigma(E_inf, E_s, eta_0, eta_inf, lam, n, i_a, epss, eps_dot)[0][-1])
        sigmas.append(sigma(lambda_inf, mu_inf, eta_0,  eta_inf, i_a, lam, n, epss, eps_dot)[0][-1])

    plt.semilogx(eps_dots, sigmas, label=r'$a = ' + str(i_a) + '$')

plt_cfig.format_and_save(f'{fig_path}/carreau-3d-as.pdf', save=save)


plt.show()