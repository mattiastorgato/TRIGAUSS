def chebyshev(N, mom, abm=None):

    if N <= 0:
        raise ValueError('N fuori range')

    # controllo dimensione rispetto a mom
    if N > mom.shape[0] / 2:
        N = mom.shape[0] // 2

    # caso in cui abm non viene fornito
    if abm is None:
        abm = np.zeros((2 * N - 1, 2))

    # controllo dimensione rispeto a abm
    if N > (abm.shape[0] + 1) / 2:
        N = int((abm.shape[0] + 1) / 2)

    # inizializza ab
    ab = np.zeros((N, 2))

    ab[0, 0] = abm[0, 0] + mom[1] / mom[0]
    ab[0, 1] = mom[0]

    # ci assicuriamo che normsq è un vettore colonna
    if N == 1:
        normsq = np.array([mom[0]])
        return ab, normsq.reshape(-1, 1)

    # Inizializziamo matrice sig
    sig = np.zeros((N + 1, 2 * N)) # righe da 0 a N, colonne da 0 a 2N-1

    # sig[0, :] rimane 0
    sig[1, :] = mom[0:2*N]

    # ciclo per il calcolo di ab
    for n in range(2, N + 1):

        for m in range(n - 3, 2 * N - n + 1):

            sig[n, m] = sig[n - 1, m + 1] + (ab[n - 2, 0] - abm[m, 0]) * sig[n - 1, m] - ab[n - 2, 1] * sig[n - 2, m] + abm[m, 1] * sig[n - 1, m - 1]


        ab[n-1, 0] = abm[n - 1, 0] + sig[n, n] / sig[n, n - 1] - sig[n - 1, n - 1] / sig[n - 1, n - 2]


        ab[n-1, 1] = sig[n, n - 1] / sig[n - 1, n - 2]

    normsq = np.zeros(N)
    for n in range(N):
        normsq[n] = sig[n + 1, n]

    return ab, normsq.reshape(-1, 1) # ci assicuriamo che normsq è un vettore colonna
