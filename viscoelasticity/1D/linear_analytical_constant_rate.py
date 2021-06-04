import numpy as np
import matplotlib.pyplot as plt


def sigma(E_inf, E_s, tau, eps, eps_dot):
    return E_inf * eps + tau * E_s * eps_dot * (1-np.exp(-eps/(eps_dot*tau)))


eps_dot = np.logspace(-4, 6, 1000)

eps = np.linspace(0, 0.01, 5)


E_inf = 2e9
E_s = 1e9
tau = 1e-3


for ep in eps[1:]:
    plt.semilogx(eps_dot, sigma(E_inf, E_s, tau, ep, eps_dot), label=r'$\varepsilon = ' + str(ep) + '$')

plt.grid()
plt.legend()
plt.xlabel(r'$\dot{\varepsilon}$')
plt.ylabel(r'$\sigma$')


plt.figure()

taus = np.logspace(-5, -1, 5)

for tau in taus:
    plt.semilogx(eps_dot, sigma(E_inf, E_s, tau, eps[-1], eps_dot), label=r'$\tau = ' + str(tau) + '\,s$')

plt.grid()
plt.legend()
plt.xlabel(r'$\dot{\varepsilon}$')
plt.ylabel(r'$\sigma$')

plt.figure()

tau = 1e-3
E_ss = np.linspace(0.8e9, 1.2e9, 5)

for E_s in E_ss:
    plt.semilogx(eps_dot, sigma(E_inf, E_s, tau, eps[-1], eps_dot), label=r'$E_s = ' + str(E_s) + '\,s$')

plt.grid()
plt.legend()
plt.xlabel(r'$\dot{\varepsilon}$')
plt.ylabel(r'$\sigma$')

plt.show()
