def gqcircsect(n, omega, r1, r2):

    # Trigonometric Gaussian formula on the arc
    tw = trigauss(n, -omega, omega)

    # Algebraic Gaussian formula on the radial segments
    m_radial = int(np.ceil((n + 2) / 2))
    ab = r_jacobi(m_radial, 0, 0)
    xw = gauss(m_radial, ab)
    xw[:, 0] = xw[:, 0] * (r2 - r1) / 2 + (r2 + r1) / 2
    xw[:, 1] = xw[:, 1] * (r2 - r1) / 2

    # Creating the polar grid
    r, theta = np.meshgrid(xw[:, 0], tw[:, 0])
    w1, w2 = np.meshgrid(xw[:, 1], tw[:, 1])

    # Nodal cartesian coordinates and weights
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    w = r * w1 * w2

    xyw = np.column_stack((x.flatten(), y.flatten(), w.flatten()))
    return xyw