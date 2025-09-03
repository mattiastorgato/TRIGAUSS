

def gauss(N, ab):

    N0 = ab.shape[0]
    if N0 < N:
        raise ValueError('l\'array di input ab è troppo corto')

    # Costruiamo la matrice di Jacobi J
    J = np.zeros((N, N))

    # Inseriamo gli alpha_k sulla diagonale principale
    for n in range(N):
        J[n, n] = ab[n, 0]

    # Inserisci le radici quadrate dei beta_k sulla sub-diagonale e super-diagonale
    for n in range(1, N):

        J[n, n - 1] = np.sqrt(ab[n, 1])

        J[n - 1, n] = J[n, n - 1]

    # Calcola autovalori (nodi) e autovettori
    autoval, autovet = np.linalg.eigh(J)

    # Ordina gli autovalori e riordina gli autovettori di conseguenza
    sorted_indices = np.argsort(autoval)
    nodes = autoval[sorted_indices]      # Ordina i nodi
    # Riordina le colonne degli autovettori in base agli indici ordinati
    autovet = autovet[:, sorted_indices]

    # Calcola i pesi
    weights = ab[0, 1] * (autovet[0, :]**2)

    xw = np.column_stack((nodes, weights))

    return xw
