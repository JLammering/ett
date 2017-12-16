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
    mu = np.logspace(-1, 30, 10000)
    alf = alpha(alpha_mz, mass_Z, mu)
    plt.fill_between(mu, noms(alf)+stds(alf), noms(alf) - stds(alf),
                     color='red', label='Ein-Sigma-Intervall')

    plt.plot(mu, noms(alf), label='Kurve')

    d = np.zeros(len(alf))
    d = noms(alf) - 1/128
    print(d)
    for i in range(len(d) - 1):
        if d[i] == 0. or d[i] * d[i + 1] < 0.:  # crossover at i
            print('Grafische Lösung:', alf[i])

    plt.xlim(0.1, 1000)
    plt.xlabel(r'$\mu \:/\: \si{\giga\electronvolt}$')
    plt.ylabel(r'$\alpha_S$')
    plt.legend(loc='best')

    plt.xscale('log')

    # in matplotlibrc leider (noch) nicht möglich
    plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
    plt.savefig('build/plotalpha.pdf')
    plt.close()


def diverg(alpha_mz, mass_Z):
    mu = unp.sqrt(unp.exp((-1)*4*np.pi/(alpha_mz*(11-(10/3))))*mass_Z**2)
    return mu


def sbeta(Nf):
    return (11 - (2/3) * Nf)


def ebeta(Nu):
    return (-16/3 * (1 + Nu/3))


def muvfunc(Nu, Nf, alpha_mz, ealpha_mz, mass_Z):
    muv = unp.exp(((ealpha_mz - alpha_mz)/(ealpha_mz*alpha_mz))*(2*np.pi/(ebeta(Nu) - sbeta(Nf))))*mass_Z
    return muv


if __name__ == '__main__':
    # z, b = np.genfromtxt('daten/magnetfeld.txt', unpack='True')
    mass_Z = ufloat(91.1876, 0.0021)
    alpha_mz = ufloat(0.1182, 0.0016)
    ealpha_mz = ufloat(127.950, 0.017)
    # ealpha_mz = ufloat(137, 0.014)
    ealpha_mz = 1/ealpha_mz
    Nf = np.array([5, 6, 10])
    Nu = np.array([2, 3, 3])
    plotalpha(alpha_mz, mass_Z)
    print('Skala divergent bei mu=', diverg(alpha_mz, mass_Z), 'GeV')
    print('muv =', muvfunc(Nu, Nf, alpha_mz, ealpha_mz, mass_Z), 'GeV')
    print(sbeta(Nf), ebeta(Nu))
