import solid_mechanics as sm
import numpy as np

# E_inf = 2e9
E_inf = 0e9
v_inf = 0.001
E_s = 1e9
v_s = 0.001
tau = 1e-3


def e_v_2_lam_mu(E, v):
    return E * v / ((1 + v) * (1 - 2 * v)), E / (2 * (1 + v))


lambda_inf, mu_inf = e_v_2_lam_mu(E_inf, v_inf)
lam_s, mu_s = e_v_2_lam_mu(E_s, v_s)

eta_d = E_s * tau

n = 0.33
eta_0 = E_s * tau
eta_inf = eta_0*1e5
a = 1

lam = 1


eps_n1 = np.array([3., 0.3, 0.3, 0., 1., 0.])
eps_dn = np.array([0., 0., 0., 0., 0., 0.])
dt = 0.1

eps_dn1, sigma = sm.solid_mechanics.increment_carraeu_visc_elastic(mu_inf,
                                                                   lambda_inf,
                                                                   lam,
                                                                   n, a,
                                                                   eta_0, eta_inf,
                                                                   mu_s, eps_n1, eps_dn, dt)

def apprx_dsig_deps(eps_n1, eps_dn, dt):
    eps_dn1, sigma = sm.solid_mechanics.increment_carraeu_visc_elastic(mu_inf,
                                                                       lambda_inf,
                                                                       lam,
                                                                       n, a,
                                                                       eta_0, eta_inf,
                                                                       mu_s, eps_n1, eps_dn, dt)

    d_eps = 1.e-7
    C = np.zeros([6, 6])
    for i in range(6):
        eps_n1_hat = eps_n1.copy()
        eps_n1_hat[i] += d_eps
        eps_dn1, sigma_hat = sm.solid_mechanics.increment_carraeu_visc_elastic(mu_inf,
                                                                               lambda_inf,
                                                                               lam,
                                                                               n, a,
                                                                               eta_0, eta_inf,
                                                                               mu_s, eps_n1_hat, eps_dn, dt)

        C[:, i] = (sigma_hat - sigma) / d_eps

    return C


proj = np.zeros([6, 6])
proj[:3, :3] = -1. / 3.
proj[range(6), range(6)] += 1
eye6 = np.eye(6, 6)

C_approx = apprx_dsig_deps(eps_n1, eps_dn, dt)
C_an_fort = sm.solid_mechanics.carraeu_visc_dsig_bar_e_deps_n1(mu_s, lam, a, n, eta_0, eta_inf, eps_n1, eps_dn, eps_dn1, dt)
a = 0
