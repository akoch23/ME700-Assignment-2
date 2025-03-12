import math
import numpy as np
import scipy
from scipy.linalg import eig
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math_utils

# Given Node positions, 3D Force, 3D Moment, b, h, E, v, J, E0, E1, A, I_y, I_z, I_rho

r = 1
E = 10000
nu = 0.3
A = np.pi * (r**2)
I_y = np.pi * ((r**4)/4)
I_z = np.pi * ((r**4)/4)
I_rho = np.pi * ((r**4)/2)
J = np.pi * ((r**4)/2)


x0, y0, z0 = 0, 0, 0
x1, y1, z1 = 25, 50, 37

P = 1
L = np.sqrt((x1 - x0)**2 + (y1-y0)**2 + (z1-z0)**2)


def define_forces_and_moments():
    F_x = -1 * P * (x1 - x0) / L
    F_y = -1 * P * (y1 - y0) / L
    F_z = -1 * P * (z1 - z0) / L

    F = [F_x, F_y, F_z]

    M_x = 0
    M_y = 0
    M_z = 0

    M = [M_x, M_y, M_z]

    return F, M


class Node:
    def __init__(self, x, y, z, constraints=None, forces=None):
        self.x = x
        self.y = y
        self.z = z
        self.constraints = constraints if constraints else [False] * 6
        self.forces = forces if forces else [0] * 6

    def __repr__(self):
        return f"Node({self.x}, {self.y}, {self.z})"    

    def coordinates(self):
        return [self.x, self.y, self.z]
    
    
class Element:
    def __init__(self, node_1, node_2, E, nu, A, I_y, I_z, I_rho, J, local_z = None):
        self.node_1 = node_1
        self.node_2 = node_2
        self.E = E
        self.nu = nu
        self.A = A
        self.I_y = I_y
        self.I_z = I_z
        self.I_rho = I_rho
        self.J = J
        self.local_z = local_z
        self.L = self.length()
    
    def __repr__(self):
        return f"Element({self.node_1}, {self.node_2})"
    
    def __eq__(self, other):
        return (self.node_1, self.node_2) == (other.node_1, other.node_2)
    
    def __hash__(self):
        return hash((self.node_1, self.node_2))
    
    def nodes(self):
        return [self.node_1, self.node_2]
    
    
    def get_forces_and_moments(self):
        # Return forces and moments at both ends (node 1 and node 2)
        Fx1, Fy1, Fz1, Mx1, My1, Mz1 = self.node_1.forces
        Fx2, Fy2, Fz2, Mx2, My2, Mz2 = self.node_2.forces
        
        print(f"\nElement {self}:")
        print(f"Node 1 Forces: {self.node_1.forces}")
        print(f"Node 2 Forces: {self.node_2.forces}")
        
        return Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2
    
    def length(self):
        x_1, y_1, z_1 = self.node_1.coordinates()
        x_2, y_2, z_2 = self.node_2.coordinates()
        
        length = math.sqrt((x_2 - x_1)**2 + (y_2 - y_1)**2 + (z_2 - z_1)**2)
        return length

    def stiffness_matrix(self):
        Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2 = self.get_forces_and_moments()
        k_g_local = local_geometric_stiffness_matrix_3D_beam(element.L, A, I_rho, Fx2, Mx2, My1, Mz1, My2, Mz2)
        k_e_local = local_elastic_stiffness_matrix_3D_beam(self.E, self.nu, self.A, self.L, self.I_y, self.I_z, self.J)
        gamma = rotation_matrix_3D(*self.node_1.coordinates(), *self.node_2.coordinates(), v_temp=self.local_z)
        transformation = transformation_matrix_3D(gamma)
        K_global = transformation.T @ k_e_local @ transformation
        return K_global

def assemble_global_elastic_stiffness_matrix(nodes, elements):
    total_dof = len(nodes) * 6
    K_e_global = np.zeros((total_dof, total_dof))

    for element in elements:
        K_e = element.stiffness_matrix()
        node_1_index = nodes.index(element.node_1)
        node_2_index = nodes.index(element.node_2)

        dof_indices = [
            node_1_index * 6, node_1_index * 6 + 1, node_1_index * 6 + 2,
            node_1_index * 6 + 3, node_1_index * 6 + 4, node_1_index * 6 + 5,
            node_2_index * 6, node_2_index * 6 + 1, node_2_index * 6 + 2,
            node_2_index * 6 + 3, node_2_index * 6 + 4, node_2_index * 6 + 5
        ]

        for i in range(12):
            for j in range(12):
                K_e_global[dof_indices[i], dof_indices[j]] += K_e[i,j]
                
    return K_e_global

def assemble_global_geometric_stiffness_matrix(nodes, elements, A, I_rho):
    total_dof = len(nodes) * 6
    K_g_global = np.zeros((total_dof, total_dof))

    forces_and_moments = {}
    
    for element in elements:
        # Extract necessary properties for each element
        node_1, node_2 = element.nodes()
        x1, y1, z1 = node_1.coordinates()
        x2, y2, z2 = node_2.coordinates()

        # Local geometric stiffness matrix for this element
        forces_and_moments_raw = element.get_forces_and_moments()  # Define or extract forces and moments
        
        Fx1 = forces_and_moments_raw[0]   # Force in X at node 1
        Fy1 = forces_and_moments_raw[1]   # Force in Y at node 1
        Fz1 = forces_and_moments_raw[2]   # Force in Z at node 1
        Mx1 = forces_and_moments_raw[3]   # Moment in X at node 1
        My1 = forces_and_moments_raw[4]   # Moment in Y at node 1
        Mz1 = forces_and_moments_raw[5]   # Moment in Z at node 1
        Fx2 = forces_and_moments_raw[6]   # Force in X at node 2
        Fy2 = forces_and_moments_raw[7]   # Force in X at node 2
        Fz2 = forces_and_moments_raw[8]   # Force in X at node 2
        Mx2 = forces_and_moments_raw[9]   # Moment in X at node 2
        My2 = forces_and_moments_raw[10]  # Moment in Y at node 2
        Mz2 = forces_and_moments_raw[11]  # Moment in Z at node 2
        
        Fx2 += Fx1
        Mx2 += Mx1
        My1 += My2
        Mz1 += Mz2
        
        forces_and_moments[element] = (Fx2, Mx2, My1, Mz1, My2, Mz2)

        
        print(f"\nForces and Moments for Element {element}:")
        print(f"Fx2: {Fx2}, Mx2: {Mx2}, My1: {My1}, Mz1: {Mz1}, My2: {My2}, Mz2: {Mz2}")

        # If all are zero, the stiffness matrix will also be zero
        if Fx2 == 0 and Mx2 == 0 and My1 == 0 and Mz1 == 0 and My2 == 0 and Mz2 == 0:
            print("Warning: Zero forces and moments detected! Geometric stiffness will be zero.")
        
        forces_and_moments[element] = (Fx2, Mx2, My1, Mz1, My2, Mz2)
        
        k_g_local = local_geometric_stiffness_matrix_3D_beam(element.L, A, I_rho, Fx2, Mx2, My1, Mz1, My2, Mz2)

        print(f"\nLocal Geometric Stiffness Matrix for Element {element}:")
        print(k_g_local)
        
        # Global DOF indices for node 1 and node 2
        node_1_index = nodes.index(node_1)
        node_2_index = nodes.index(node_2)

        # Global indices for DOFs associated with the two nodes
        dof_indices_node_1 = [node_1_index * 6 + i for i in range(6)]
        dof_indices_node_2 = [node_2_index * 6 + i for i in range(6)]

        # Assemble the global geometric stiffness matrix by adding the contributions from each element
        for i in range(12):
            for j in range(12):
                # Mapping the local stiffness matrix to the global matrix
                global_i = dof_indices_node_1[i] if i < 6 else dof_indices_node_2[i - 6]
                global_j = dof_indices_node_1[j] if j < 6 else dof_indices_node_2[j - 6]
                K_g_global[global_i, global_j] += k_g_local[i, j]
                
    print("\nGlobal Geometric Stiffness Matrix (K_g):")
    print(K_g_global)
    print("Sum of K_g elements:", np.sum(K_g))

    return K_g_global

def apply_boundary_conditions(K, F, nodes):
    constrained_DOFs = []
    for i, node in enumerate(nodes):
        for dof in range(6):
            if node.constraints[dof]:
                constrained_DOFs.append(i * 6 + dof)
                
    constrained_DOFs.sort()
    K_mod = np.delete(K, constrained_DOFs, axis = 0)
    K_mod = np.delete(K_mod, constrained_DOFs, axis=1)
    F_mod = np.delete(F, constrained_DOFs, axis=0)

    return K_mod, F_mod, constrained_DOFs

def solve_displacements(K_mod, F_mod):
    return np.linalg.solve(K_mod, F_mod)

def compute_reaction_forces(K, U_full, constrained_DOFs):
    F_internal = K @ U_full
    return F_internal


def node_coordinates():
    node_num = 7
    node_coords = []
    
    for i in range(node_num):        
        print(f"Enter coordinates for Node {i}:")
        
        x = float(input(f"  x-coordinate for Node {i}: "))
        y = float(input(f"  y-coordinate for Node {i}: "))
        z = float(input(f"  z-coordinate for Node {i}: "))

        print("Select node constraint type: \n 1 - Fixed Node \n 2 - Pinned Node \n 3 - Free Node \n 4 - Fixed Node (X only) \n 5 - Fixed Node (Y only) \n 6 - Fixed Node (Z only) \n 7 - Fixed Node (X-rotation only) \n 8 - Fixed Node (Y-rotation only) \n 9 - Fixed Node (Z-rotation only)")
        node_type = int(input(f" Enter the number corresponding to constraint type for Node {i}: "))

        if node_type == 1:
            constraints = [True, True, True, True, True, True]
        elif node_type == 2:
            constraints = [True, True, True, False, False, False]
        elif node_type == 3:
            constraints = [False, False, False, False, False, False]
        elif node_type == 4:
            constraints = [True, False, False, False, False, False]
        elif node_type == 5:
            constraints = [False, True, False, False, False, False]
        elif node_type == 6:
            constraints = [False, False, True, False, False, False]
        elif node_type == 7:
            constraints = [False, False, False, True, False, False]
        elif node_type == 8:
            constraints = [False, False, False, False, True, False]
        elif node_type == 9:
            constraints = [False, False, False, False, False, True]
        
        new_node = Node(x, y, z, constraints=constraints)
        node_coords.append(new_node)
    
    return node_coords
    

def apply_forces_and_moments(node_pos, F, M):

    applied_forces = {i: np.zeros(6) for i in range(len(node_pos))}

    force_nodes = input("Enter the node index where the force vector (F) should be applied (space-separated): ")
    force_nodes = list(map(int, force_nodes.split())) if force_nodes.strip() else []

    moment_nodes = input("Enter the node index where the moment vector (M) should be applied (space-separated): ")
    moment_nodes = list(map(int, moment_nodes.split())) if moment_nodes.strip() else []

    for i in force_nodes:
        if 0 <= i < len(node_pos):
            applied_forces[i][:3] = F[:]
            node_pos[i].forces[:3] = F[:]

    for i in moment_nodes:
        if 0 <= i < len(node_pos):
            applied_forces[i][3:] = M[:]
            node_pos[i].forces[3:] = M[:]
            
    for node, force in applied_forces.items():
        print(f"Node {node} - Applied Force: {force[:3]}, Applied Moment: {force[3:]}")


    return applied_forces



F, M = define_forces_and_moments()

node_pos = node_coordinates()

applied_loads = apply_forces_and_moments(node_pos, F, M)


print("Defined Nodal Positions:")
for i, node in enumerate(node_pos, start=1):
    print(f"Node {i}: {node.coordinates()}")

for node_idx, load in applied_loads.items():
    print(f"Node {node_idx}: Applied Forces = {load[:3]}, Applied Moments = {load[3:]}")
    
def define_elements(nodes):
    
    element_num = 6
    element_connect = []
    
    for i in range(element_num):
        print(f"Enter the nodes that serve as the end points for Element {i}:")
        
        node_1_index = int(input(f" Enter the index of Node 1 for Element {i}: "))
        node_2_index = int(input(f" Enter the index of Node 2 for Element {i}: "))

        local_z_input = input(f"Do you want to define a local_z vector for Element {i}? (y/n): ").strip().lower()

        if local_z_input == 'y':
            print(f"Enter the local_z vector for Element {i} (as x, y, z components):")
            x = float(input(f"  x-component of local_z for Element {i}: "))
            y = float(input(f"  y-component of local_z for Element {i}: "))
            z = float(input(f"  z-component of local_z for Element {i}: "))
            local_z = np.array([x, y, z])

        else:
            local_z = None
        
        element = Element(nodes[node_1_index], nodes[node_2_index], E, nu, A, I_y, I_z, I_rho, J, local_z)
        element_connect.append(element)
        
    return element_connect

elements = define_elements(node_pos)

K_e_global = assemble_global_elastic_stiffness_matrix(node_pos, elements)
print(K_e_global)
K_g_global = assemble_global_geometric_stiffness_matrix(node_pos, elements, A, I_rho)

print("\nFinal Global Geometric Stiffness Matrix:")
print(K_g_global)
print("Sum of K_g elements:", np.sum(K_g))

print("Defined Elements and Connected Nodes:")
for i, element in enumerate(elements, start=1):
    node_1, node_2 = element.nodes()
    length = element.length()
    print(f"Element {i}: connects Node({node_1.coordinates()}) and Node({node_2.coordinates()}) with length {length: .1f}")

print("Global Elastic Stiffness Matrix: ", K_e_global)


F_global = np.zeros(len(node_pos) * 6)

for node_idx, load in applied_loads.items():
    dof_index = node_idx * 6
    F_global[dof_index:dof_index+6] = load  # Assign applied forces and moments

K_mod, F_mod, constrained_DOFs = apply_boundary_conditions(K_e_global, F_global, node_pos)

# Solve for displacements
U_mod = solve_displacements(K_mod, F_mod)

# Reconstruct full displacement vector (including constrained DOFs)
U_full = np.zeros(len(node_pos) * 6)
free_dof_indices = [i for i in range(len(node_pos) * 6) if i not in constrained_DOFs]

for i, dof in enumerate(free_dof_indices):
    U_full[dof] = U_mod[i]

    
def plot_deformed_shape(nodes, elements, U_full, scale=1.0):
    """
    Plot the undeformed and deformed shape of the structure in 3D.
    
    Parameters:
        nodes: List of Node objects.
        elements: List of Element objects.
        U_full: Displacement vector (including constrained DOFs).
        scale: Scaling factor for visualization.
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for element in elements:
        node_1, node_2 = element.nodes()

        # Get original coordinates
        x = [node_1.x, node_2.x]
        y = [node_1.y, node_2.y]
        z = [node_1.z, node_2.z]
        
        # Plot undeformed structure (original beam position)
        ax.plot(x, y, z, 'bo-', label='Undeformed' if element == elements[0] else "")

        # Get displacement indices
        idx1 = nodes.index(node_1) * 6
        idx2 = nodes.index(node_2) * 6

        # Compute deformed coordinates
        xd = [x[0] + scale * U_full[idx1], x[1] + scale * U_full[idx2]]
        yd = [y[0] + scale * U_full[idx1+1], y[1] + scale * U_full[idx2+1]]
        zd = [z[0] + scale * U_full[idx1+2], z[1] + scale * U_full[idx2+2]]

        # Plot deformed structure (scaled displacement)
        ax.plot(xd, yd, zd, 'r--', label='Deformed' if element == elements[0] else "")

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Deformed Shape of the Structure')
    ax.legend()
    plt.show()

def plot_internal_forces(elements, internal_forces):
    """
    Plot internal forces and moments for each element in local coordinates.
    
    Parameters:
        elements: List of Element objects.
        internal_forces: Dictionary of internal forces/moments for each element.
    """
    for element in elements:
        x = [0, element.L]  # Local element length
        forces = internal_forces[element]

        fig, ax = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle(f'Internal Forces for Element {element}', fontsize=16)

        # Axial Force
        ax[0, 0].plot(x, forces[[0, 6]], 'r', label='Axial Force')
        ax[0, 0].set_title('Axial Force')
        ax[0, 0].grid()

        # Shear Force Y
        ax[0, 1].plot(x, forces[[1, 7]], 'g', label='Shear Force Y')
        ax[0, 1].set_title('Shear Force Y')
        ax[0, 1].grid()

        # Shear Force Z
        ax[1, 0].plot(x, forces[[2, 8]], 'b', label='Shear Force Z')
        ax[1, 0].set_title('Shear Force Z')
        ax[1, 0].grid()

        # Bending Moment Z
        ax[1, 1].plot(x, forces[[5, 11]], 'm', label='Bending Moment Z')
        ax[1, 1].set_title('Bending Moment Z')
        ax[1, 1].grid()

        plt.show()
    
# Compute reaction forces
reaction_forces = compute_reaction_forces(K_e_global, U_full, constrained_DOFs)

plot_deformed_shape(node_pos, elements, U_full, scale=1000)

internal_forces = {}
for element in elements:
    node_1, node_2 = element.nodes()
    idx1 = node_pos.index(node_1) * 6
    idx2 = node_pos.index(node_2) * 6
    dof_indices = np.array([
        idx1, idx1+1, idx1+2, idx1+3, idx1+4, idx1+5,
        idx2, idx2+1, idx2+2, idx2+3, idx2+4, idx2+5
    ])
    
    # Get displacements for the element
    d_elem = U_full[dof_indices]
    
    # Compute local internal forces using element stiffness
    k_local = element.stiffness_matrix()
    internal_forces[element] = k_local @ d_elem

plot_internal_forces(elements, internal_forces)

# Display results
print("\nNodal Displacements:")
for i in range(len(node_pos)):
    print(f"Node {i}: Displacements = {U_full[i*6:i*6+3]}, Rotations = {U_full[i*6+3:i*6+6]}")

# Compute total force/moment at each node (includes both applied and reaction forces)
total_forces = np.zeros(len(node_pos) * 6)

for i in range(len(node_pos) * 6):
    if i in constrained_DOFs:
        total_forces[i] = reaction_forces[i]  # Assign reaction force at constrained DOFs
    else:
        total_forces[i] = F_global[i]  # Assign applied force at non-constrained DOFs

# Display results in a node-based format
print("\nNodal Forces and Moments:")
for i in range(len(node_pos)):
    print(f"Node {i}: Forces = {total_forces[i*6:i*6+3]}, Moments = {total_forces[i*6+3:i*6+6]}")


def compute_elastic_critical_load_factor(K_e, K_g):
    
    print("\nElastic Stiffness Matrix (K_e):")
    print(K_e)
    print("Sum of K_e elements:", np.sum(K_e))
    
    print("\nGeometric Stiffness Matrix (K_g):")
    print(K_g)
    print("Sum of K_g elements:", np.sum(K_g))
    
    # Check for singularity
    if np.linalg.cond(K_g) > 1e12:  # Arbitrary large number indicating near singularity
        print("Warning: K_g is nearly singular!")
    
    # Solve the generalized eigenvalue problem K_e * u = λ * K_g * u
    eigenvalues, eigenvectors = eig(K_e, K_g)
    print(f"Eigenvalues: {eigenvalues}")
    
    # Sort eigenvalues in ascending order
    eigenvalues = np.real(eigenvalues)
    eigenvalues.sort()
    
    # The smallest eigenvalue is the Elastic Critical Load Factor (λ)
    elastic_critical_load_factor = eigenvalues[0]
    return elastic_critical_load_factor

# Assuming you have the global stiffness matrices K_e and K_g
K_e = K_e_global  # The global elastic stiffness matrix
K_g = K_g_global

    
elastic_critical_load_factor = compute_elastic_critical_load_factor(K_e, K_g)

print(f"Elastic Critical Load Factor (λ): {elastic_critical_load_factor}")
