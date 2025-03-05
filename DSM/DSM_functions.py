import numpy as np
# Given Node positions, 3D Force, 3D Moment, b, h, E, v, J, E0, E1, A, I_y, I_z, I_rho


def local_element_stiffness_matrix(E: float, nu: float, A: float, L: float, Iy: float, Iz: float, J: float) -> np.ndarray:
    k = np.zeros((12, 12))
    
    k[0, 0] = A/L
    k[0, 6] = -A/L
    k[6, 0] = -A/L
    k[6, 6] = A/L

    k[3, 3] = (12 * I_y)/L**3
    k_e[3, 8] = -(12 * I_y)/L**3  
    k_e[8, 3] = -(12 * I_y)/L**3 
    k_e[8, 8] = (12 * I_y)/L**3 

    k[1, 1] = (12 * I_z)/L**3
    k[1, 7] = -(12 * I_z)/L**3
    k[7, 1] = -(12 * I_z)/L**3
    k[7, 7] = (12 * I_z)/(L**3)
    k[1, 5] = (6 * I_z)/ L**2
    k[5, 1] = (6 * I_z)/ L**2
    k[1, 11] = (6 * I_z)/L**2
    k[11, 1] = (6 * I_z)/L**2
    k[5, 7] = -(6 * I_z)/L**2
    k[7, 5] = -(6 * I_z)/L**2
    k[7, 11] = -(6 * I_z)/L**2
    k[11, 7] = -(6 * I_z)/L**2
    k[5, 5] = (4 * I_z)/L
    k[11, 11] = (4 * I_z)/L
    k_e[5, 11] = (2 * I_z)/L
    k_e[11, 5] = (2 * I_z)/L

    k[2, 2] = (12 * I_y)/L**3
    k[2, 8] = -(12 * I_y)/L**3
    k[8, 2] = -(12 * I_y)/L**3
    k[8, 8] = (12 * I_y)/L**3
    k[2, 4] = -(6 * I_y)/L**2
    k[4, 2] = -(6 * I_y)/L**2
    k[2, 10] = -(6 * I_y)/L**2
    k[10, 2] = -(6 * I_y)/L**2
    k[4, 8] = (6 * I_y)/L**2
    k[8, 4] = (6 * I_y)/L**2
    k[8, 10] = (6 * I_y)/L**2
    k[10, 8] = (6 * I_y)/L**2
    k[4, 4] = (4 * I_y)/L
    k[10, 10] = (4 * I_y)/L
    k[4, 10] = (2 * I_y)/L
    k[10, 4] = (2 * I_y)/L
    return k


def rotation_matrix(x1: float, y1: float, z1: float, x2: float, y2: float, z2: float, z: np.ndarray = None):
    L = np.sqrt((x2 - x1) ** 2.0 + (y2 - y1) ** 2.0 + (z2 - z1) ** 2.0)
    lxp = (x2 - x1) / L
    mxp = (y2 - y1) / L
    nxp = (z2 - z1) / L
    local_x = np.asarray([lxp, mxp, nxp])

    v_temp = z
    # choose a vector to make orthogonal relative to the y axis if one is not given
    if v_temp is None:
        if np.isclose(lxp, 0.0) and np.isclose(mxp, 0.0):
            v_temp = np.array([0, 1.0, 0.0])
        else:
            # otherwise use the global z axis
            v_temp = np.array([0, 0, 1.0])
    else:
        # check to make sure that given v_temp is a unit vector
        check_unit_vector(v_temp)
        # check to make sure that given v_temp is not parallel to the local x axis
        check_parallel(local_x, v_temp)
    
    local_y = np.cross(v_temp, local_x)
    local_y = local_y / np.linalg.norm(local_y)
    local_z = np.cross(local_x, local_y)
    local_z = local_z / np.linalg.norm(local_z)
    gamma = np.vstack((local_x, local_y, local_z))
    
    return gamma


def transformation_matrix(gamma: np.ndarray) -> np.ndarray:
    Gamma = np.zeros((12, 12))
    Gamma[0:3, 0:3] = gamma
    Gamma[3:6, 3:6] = gamma
    Gamma[6:9, 6:9] = gamma
    Gamma[9:12, 9:12] = gamma
    return Gamma


def local_geometric_stiffness_matrix(L, A, I_rho, Fx2, Mx2, My1, Mz1, My2, Mz2):
    k_g = np.zeros((12, 12))
    # upper triangle off diagonal terms
    k_g[0, 6] = -Fx2 / L
    k_g[1, 3] = My1 / L
    k_g[1, 4] = Mx2 / L
    k_g[1, 5] = Fx2 / 10.0
    k_g[1, 7] = -6.0 * Fx2 / (5.0 * L)
    k_g[1, 9] = My2 / L
    k_g[1, 10] = -Mx2 / L
    k_g[1, 11] = Fx2 / 10.0
    k_g[2, 3] = Mz1 / L
    k_g[2, 4] = -Fx2 / 10.0
    k_g[2, 5] = Mx2 / L
    k_g[2, 8] = -6.0 * Fx2 / (5.0 * L)
    k_g[2, 9] = Mz2 / L
    k_g[2, 10] = -Fx2 / 10.0
    k_g[2, 11] = -Mx2 / L
    k_g[3, 4] = -1.0 * (2.0 * Mz1 - Mz2) / 6.0
    k_g[3, 5] = (2.0 * My1 - My2) / 6.0
    k_g[3, 7] = -My1 / L
    k_g[3, 8] = -Mz1 / L
    k_g[3, 9] = -Fx2 * I_rho / (A * L)
    k_g[3, 10] = -1.0 * (Mz1 + Mz2) / 6.0
    k_g[3, 11] = (My1 + My2) / 6.0
    k_g[4, 7] = -Mx2 / L
    k_g[4, 8] = Fx2 / 10.0
    k_g[4, 9] = -1.0 * (Mz1 + Mz2) / 6.0
    k_g[4, 10] = -Fx2 * L / 30.0
    k_g[4, 11] = Mx2 / 2.0
    k_g[5, 7] = -Fx2 / 10.0
    k_g[5, 8] = -Mx2 / L
    k_g[5, 9] = (My1 + My2) / 6.0
    k_g[5, 10] = -Mx2 / 2.0
    k_g[5, 11] = -Fx2 * L / 30.0
    k_g[7, 9] = -My2 / L
    k_g[7, 10] = Mx2 / L
    k_g[7, 11] = -Fx2 / 10.0
    k_g[8, 9] = -Mz2 / L
    k_g[8, 10] = Fx2 / 10.0
    k_g[8, 11] = Mx2 / L
    k_g[9, 10] = (Mz1 - 2.0 * Mz2) / 6.0
    k_g[9, 11] = -1.0 * (My1 - 2.0 * My2) / 6.0
    # add in the symmetric lower triangle
    k_g = k_g + k_g.transpose()
    # add diagonal terms
    k_g[0, 0] = Fx2 / L
    k_g[1, 1] = 6.0 * Fx2 / (5.0 * L)
    k_g[2, 2] = 6.0 * Fx2 / (5.0 * L)
    k_g[3, 3] = Fx2 * I_rho / (A * L)
    k_g[4, 4] = 2.0 * Fx2 * L / 15.0
    k_g[5, 5] = 2.0 * Fx2 * L / 15.0
    k_g[6, 6] = Fx2 / L
    k_g[7, 7] = 6.0 * Fx2 / (5.0 * L)
    k_g[8, 8] = 6.0 * Fx2 / (5.0 * L)
    k_g[9, 9] = Fx2 * I_rho / (A * L)
    k_g[10, 10] = 2.0 * Fx2 * L / 15.0
    k_g[11, 11] = 2.0 * Fx2 * L / 15.0
    return k_g


def local_geometric_stiffness_matrix_3D_beam_without_interaction_terms(L, A, I_rho, Fx2):
    k_g = np.zeros((12, 12))
    # upper triangle off diagonal terms
    k_g[0, 6] = -Fx2 / L
    k_g[1, 5] = Fx2 / 10.0
    k_g[1, 7] = -6.0 * Fx2 / (5.0 * L)
    k_g[1, 11] = Fx2 / 10.0
    k_g[2, 4] = -Fx2 / 10.0
    k_g[2, 8] = -6.0 * Fx2 / (5.0 * L)
    k_g[2, 10] = -Fx2 / 10.0
    k_g[3, 9] = -Fx2 * I_rho / (A * L)
    k_g[4, 8] = Fx2 / 10.0
    k_g[4, 10] = -Fx2 * L / 30.0
    k_g[5, 7] = -Fx2 / 10
    k_g[5, 11] = -Fx2 * L / 30.0
    k_g[7, 11] = -Fx2 / 10.0
    k_g[8, 10] = Fx2 / 10.0
    # add in the symmetric lower triangle
    k_g = k_g + k_g.transpose()
    # add diagonal terms
    k_g[0, 0] = Fx2 / L
    k_g[1, 1] = 6.0 * Fx2 / (5.0 * L)
    k_g[2, 2] = 6.0 * Fx2 / (5.0 * L)
    k_g[3, 3] = Fx2 * I_rho / (A * L)
    k_g[4, 4] = 2.0 * Fx2 * L / 15.0
    k_g[5, 5] = 2.0 * Fx2 * L / 15.0
    k_g[6, 6] = Fx2 / L
    k_g[7, 7] = 6.0 * Fx2 / (5.0 * L)
    k_g[8, 8] = 6.0 * Fx2 / (5.0 * L)
    k_g[9, 9] = Fx2 * I_rho / (A * L)
    k_g[10, 10] = 2.0 * Fx2 * L / 15.0
    k_g[11, 11] = 2.0 * Fx2 * L / 15.0
    return k_g
