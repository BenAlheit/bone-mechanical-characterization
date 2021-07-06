import numpy as np
import matplotlib.pyplot as plt
import config.plotting_config

E_inf = 2e9
E_s = 1e9
tau = 1e-3

eps_dot = np.logspace(-4, 6, 1000)
eps = np.linspace(0.001, 0.01, 5)

epss = np.linspace(0, 0.01, 500)
eps_dots = np.logspace(-4, 6, 10)
E_ss = np.linspace(0.8e9, 1.2e9, 5)
taus = np.logspace(-5, -1, 5)

fig_path = '../../written-work/viscoelasticity/scratch/figures/'


def sigma(E_inf, E_s, tau, eps, eps_dot):
    return E_inf * eps + tau * E_s * eps_dot * (1-np.exp(-eps/(eps_dot*tau)))


def format_and_save(path, x_label=r'$\dot{\varepsilon}$'):
    plt.grid()
    plt.legend()
    plt.xlabel(x_label)
    plt.ylabel(r'$\sigma$')
    plt.tight_layout()
    plt.savefig(path)

def main():

    for ep in eps:
        plt.semilogx(eps_dot, sigma(E_inf, E_s, tau, ep, eps_dot), label=r'$\varepsilon = ' + str(ep) + '$')

    format_and_save(f'{fig_path}e-dot-stress.pdf')
    plt.figure()

    for i_tau in taus:
        plt.semilogx(eps_dot, sigma(E_inf, E_s, i_tau, eps[-1], eps_dot), label=r'$\tau = ' + str(i_tau) + '\,s$')

    format_and_save(f'{fig_path}e-dot-stress-tau.pdf')
    plt.figure()

    for i_E_s in E_ss:
        plt.semilogx(eps_dot, sigma(E_inf, i_E_s, tau, eps[-1], eps_dot), label=r'$E_s = ' + str(i_E_s) + '\,s$')

    format_and_save(f'{fig_path}e-dot-stress-E.pdf')
    plt.figure()

    for i_eps_dot in eps_dots:
        plt.plot(epss, sigma(E_inf, E_s, tau, epss, i_eps_dot), label=r'$\dot{\varepsilon} = ' + str(i_eps_dot) + '\,/s$')

    format_and_save(f'{fig_path}e-stress.pdf', r'$\varepsilon$')

    plt.show()


if __name__ == '__main__':
    main()
