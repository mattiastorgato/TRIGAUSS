import numpy as np
from scipy.integrate import quad

def trigauss(n, alpha, beta):

    omega = (beta - alpha) / 2


    # ricorsione dei momenti di Chebyshev modificati
    z = np.zeros(n + 2)
    z[0] = 2 * omega

    integral_result, _ = quad(lambda t: np.cos(2 * n * np.arccos(np.sin(t / 2) / np.sin(omega / 2))),-omega, omega, limit=5000)
    z[n] = integral_result

    #vettore ausiliare
    temp = np.arange(2, 2 * n, 2)
    #vettori di ricorsione dei momenti
    dl = 1 / 4 - 1 / (4 * (temp - 1))
    dc = 1 / 2 - 1 / np.sin(omega / 2)**2 - 1 / (2 * (temp**2 - 1))
    du = 1 / 4 + 1 / (4 * (temp + 1))

    d =( 4 * np.cos(omega / 2) / np.sin(omega / 2)) / (temp**2 - 1)

    #caso in cui n non è abbastaza lungo
    if n > 1:
        d[n - 2] = d[n - 2] - du[n - 2] * z[n]
        z[1:n] = tridisolve(dl[1:n-1], dc[0:n-1], du[0:n-2], d[0:n-1])

    #vettore dei momenti
    mom = np.zeros(2 * n + 2)
    mom[0:2*n+2:2] = z[0:n+1]

    # normalizzazione dei momenti
    k_vals = np.arange(2, len(mom))+1;
    mom[2:] = np.exp((2 - k_vals) * np.log(2)) * mom[2:]

    # coefficienti di ricorsione dei polinomi ortogonali
    abm = np.zeros((2 * n + 1, 2))
    abm[:, 1] = 0.25
    abm[0, 1] = np.pi
    abm[1, 1] = 0.5

    ab, normsq = chebyshev(n + 1, mom, abm)

    # formula gaussiana per i pesi e nodi
    xw = gauss(n + 1, ab)

    # angoli e pesi per la formula trigonometrica
    tw = np.zeros((n + 1, 2))
    tw[:, 0] = 2 * np.arcsin(np.sin(omega / 2) * xw[:, 0]) + (beta + alpha) / 2
    tw[:, 1] = xw[:, 1]

    return tw


def tridisolve(a, b, c, d):

    x = np.copy(d)
    n = len(x)

    for j in range(n - 1):
        mu = a[j] / b[j]
        b[j + 1] = b[j + 1] - mu * c[j]
        x[j + 1] = x[j + 1] - mu * x[j]

    x[n - 1] = x[n - 1] / b[n - 1]
    for j in range(n - 2, -1, -1):
        x[j] = (x[j] - c[j] * x[j + 1]) / b[j]
    return x

