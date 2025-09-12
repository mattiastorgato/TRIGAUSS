def gqcircsegm(n, omega, r):
    tw = trigauss(n + 2, -omega, omega)

    m_r = int(np.ceil((n + 1) / 2))
    ab = r_jacobi(m_r, 0, 0)
    xw = gauss(m_r, ab)

    T, Theta = np.meshgrid(xw[:, 0], tw[1:m_r, 0], indexing='ij')
    W1, W2 = np.meshgrid(xw[:, 1], tw[1:m_r, 1], indexing='ij')

    s = np.sin(Theta).ravel()
    X = (r * np.cos(Theta)).ravel()
    Y = (r * T.ravel() * s).ravel()
    W = (r**2 * (s**2) * W1.ravel() * W2.ravel())

    xyw = np.column_stack([X, Y, W])
    return xyw