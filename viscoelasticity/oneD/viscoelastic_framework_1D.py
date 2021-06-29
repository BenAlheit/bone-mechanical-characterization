from abc import ABC, abstractmethod
import numpy as np
import scipy.optimize as op
from scipy.interpolate import lagrange


class VEConstitutiveModel(ABC):

    @abstractmethod
    def sig_d(self, eps_d_dot):
        raise NotImplementedError

    @abstractmethod
    def dsig_d_deps_d_dot(self, eps_d_dot):
        raise NotImplementedError

    @abstractmethod
    def sig_e(self, eps_e):
        raise NotImplementedError

    @abstractmethod
    def dsig_e_deps_e(self, eps_e):
        raise NotImplementedError


class LinearVEMaterial(VEConstitutiveModel):

    def __init__(self, Es, tau):
        self._eta = Es * tau
        self._Es = Es

    def sig_d(self, eps_d_dot):
        return self._eta * eps_d_dot

    def dsig_d_deps_d_dot(self, eps_d_dot):
        return self._eta

    def sig_e(self, eps_e):
        return self._Es * eps_e

    def dsig_e_deps_e(self, eps_e):
        return self._Es


class LinearEOdWVMaterial(VEConstitutiveModel):

    def __init__(self, Es, K, n):
        self._K = K
        self._n = n
        self._Es = Es

    def sig_d(self, eps_d_dot):
        return self._K * np.sign(eps_d_dot) * np.abs(eps_d_dot) ** self._n

    def dsig_d_deps_d_dot(self, eps_d_dot):
        return self._eta

    def sig_e(self, eps_e):
        return self._Es * eps_e

    def dsig_e_deps_e(self, eps_e):
        return self._Es

class LinearEMOdWVMaterial(VEConstitutiveModel):

    def __init__(self, Es, K, n):
        self._K = K
        self._n = n
        self._Es = Es

    def sig_d(self, eps_d_dot):
        return self._K * (1 + np.abs(eps_d_dot)) ** (self._n-1) * eps_d_dot

    def dsig_d_deps_d_dot(self, eps_d_dot):
        return self._eta

    def sig_e(self, eps_e):
        return self._Es * eps_e

    def dsig_e_deps_e(self, eps_e):
        return self._Es


class LinearECarreuVMaterial(VEConstitutiveModel):

    def __init__(self, Es, eta_0, eta_inf, lam, n, a):
        self._eta_0 = eta_0
        self._eta_inf = eta_inf
        self._lam = lam
        self._a = a
        self._n = n
        self._Es = Es

    def sig_d(self, eps_d_dot):
        return eps_d_dot * (self._eta_inf +
                            (self._eta_0 - self._eta_inf) *
                            (1 + (self._lam * np.abs(eps_d_dot)) ** self._a) ** ((self._n - 1) / self._a))

    def dsig_d_deps_d_dot(self, eps_d_dot):
        # return self._eta
        pass

    def sig_e(self, eps_e):
        return self._Es * eps_e

    def dsig_e_deps_e(self, eps_e):
        return self._Es


class TemporalApproximation(ABC):

    def __init__(self, n):
        self._n = n
        self._tns = np.zeros(n)
        # self._tns = np.linspace(-1, 0, n)
        self._ans = np.zeros(n)
        self._n_t_updates = 0

    # def a(self, t):
    #     raise NotImplementedError
    #
    # def a_dot(self, t):
    #     raise NotImplementedError

    @abstractmethod
    def a_dot_n1(self, a_n1):
        raise NotImplementedError

    @abstractmethod
    def da_dot_dan1(self, a_n1):
        raise NotImplementedError

    def update_tns(self, tn1):
        self._tns[1:] = self._tns[:-1]
        self._tns[0] = tn1
        self._n_t_updates += 1
        if self._n_t_updates > self._n:
            self._n_t_updates = self._n

    def update_ans(self, an1):
        self._ans[1:] = self._ans[:-1]
        self._ans[0] = an1

    def previous_a(self):
        return self._ans[1]

    def a(self):
        return self._ans[0]


class BackwardEuler(TemporalApproximation):

    def __init__(self):
        super().__init__(2)

    def a_dot_n1(self, a_n1):
        return (a_n1 - self._ans[1]) / (self._tns[0] - self._tns[1])

    def da_dot_dan1(self, a_n1):
        return 1 / (self._tns[0] - self._tns[1])


class Lagrange(TemporalApproximation):

    def __init__(self, n):
        super().__init__(n)

    def a_dot_n1(self, a_n1):
        poly = lagrange(self._tns[:self._n_t_updates], [a_n1] + self._ans.tolist()[:self._n_t_updates])
        return poly.deriv(1)(self._tns[0])

    def da_dot_dan1(self, a_n1):
        pass
        # return 1 / (self._tns[0] - self._tns[1])


class VE1D:

    def __init__(self,
                 constitutive_model: VEConstitutiveModel,
                 temporal_approximation: TemporalApproximation):
        self.cm = constitutive_model
        self.ta = temporal_approximation
        self._sig_e = 0

    def increment_model(self, t_n1, eps_n1):
        self.ta.update_tns(t_n1)

        def r(eps_d_n1):
            return self.cm.sig_d(self.ta.a_dot_n1(eps_d_n1)) - self.cm.sig_e(eps_n1 - eps_d_n1)

        eps_d_n1 = op.root(r, self.ta.previous_a()).x[0]
        self.ta.update_ans(eps_d_n1)

        self._sig_e = self.cm.sig_e(eps_n1 - eps_d_n1)

    @property
    def sig_e(self):
        return self._sig_e
