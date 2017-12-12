import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit, fmin
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from uncertainties import ufloat


def beta():
    N_f = 5
    return(11 - 2/3 * N_f)


def alpha(alpha_mz, mass_Z, mu):
    alf = alpha_mz/(1 + alpha_mz/(4*np.pi) * beta()
                    * unp.log((mu**2)/(mass_Z**2)))

    return(alf)


def plotalpha(alpha_mz, mass_Z):
    mu = np.logspace(-1, 3, 10000)
    alf = alpha(alpha_mz, mass_Z, mu)
    plt.fill_between(mu, noms(alf)+stds(alf), noms(alf) - stds(alf),
                     color='red', label='Ein-Sigma-Intervall')

    plt.plot(mu, noms(alf), label='Kurve')

    plt.xlim(0.1, 1000)
    plt.xlabel(r'$\mu \:/\: \si{\giga\electronvolt}$')
    plt.ylabel(r'$\alpha_S$')
    plt.legend(loc='best')

    plt.xscale('log')

    # in matplotlibrc leider (noch) nicht möglich
    plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
    plt.savefig('build/plotalpha.pdf')
    plt.close()


if __name__ == '__main__':
    # z, b = np.genfromtxt('daten/magnetfeld.txt', unpack='True')
    mass_Z = ufloat(91.1876, 0.0021)
    alpha_mz = ufloat(0.1181, 0.0011)
    print(alpha_mz/(1 + alpha_mz/(4*np.pi) * beta() *
                    unp.log((0.1**2)/(mass_Z**2))))
    plotalpha(alpha_mz, mass_Z)
