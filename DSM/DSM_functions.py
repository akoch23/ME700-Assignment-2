import math
import numpy as np
%run DSM_math_util.ipynb

# Given Node positions, 3D Force, 3D Moment, b, h, E, v, J, E0, E1, A, I_y, I_z, I_rho

b = float(input("Enter the base value for each element: "))
h = float(input("Enter the height value for each element: "))
A = b*h
E = float(input("Enter Young's Modulus (E) for each element: ")) 
nu = float(input("Enter Poisson's ratio (v) for each element: "))
I_y = h*(b**3/12)
I_z = b*(h**3/12)
I_rho = b*(h/12)*(b**2 + h**2)
J = float(input("Enter torsional constant (J) for each element: "))

def define_forces_and_moments():
    F_x = float(input("Normal force acting in x-direction: ")) 
    F_y = float(input("Normal force acting in y-direction: ")) 
    F_z = float(input("Normal force acting in z-direction: ")) 

    F = [F_x, F_y, F_z]

    M_x = float(input("Moment force acting in x-direction: ")) 
    M_y = float(input("Moment force acting in y-direction: ")) 
    M_z = float(input("Moment force acting in z-direction: ")) 

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
    def __init__(self, node_1, node_2, E, nu, A, I_y, I_z, J, local_z = None):
        self.node_1 = node_1
        self.node_2 = node_2
        self.E = E
        self.nu = nu
        self.A = A
        self.I_y = I_y
        self.I_z = I_z
        # self.I_rho = I_rho
        self.J = J
        self.local_z = local_z
        self.L = self.length()
    
    def __repr__(self):
        return f"Element({self.node_1}, {self.node_2})"
    
    def nodes(self):
        return [self.node_1, self.node_2]
    
    def length(self):
        x_1, y_1, z_1 = self.node_1.coordinates()
        x_2, y_2, z_2 = self.node_2.coordinates()
        
        length = math.sqrt((x_2 - x_1)**2 + (y_2 - y_1)**2 + (z_2 - z_1)**2)
        return length

    def stiffness_matrix(self):
        K_local = local_elastic_stiffness_matrix_3D_beam(self.E, self.nu, self.A, self.L, self.I_y, self.I_z, self.J)
        gamma = rotation_matrix_3D(*self.node_1.coordinates(), *self.node_2.coordinates(), v_temp=self.local_z)
        transformation = transformation_matrix_3D(gamma)
        K_global = transformation.T @ K_local @ transformation
        return K_global

def assemble_global_stiffness_matrix(nodes, elements):
    total_dof = len(nodes) * 6
    K_global = np.zeros((total_dof, total_dof))

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
                K_global[dof_indices[i], dof_indices[j]] += K_e[i,j]
    return K_global

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
    return K @ U_full


def node_coordinates():
    node_num = int(input("Enter the number of nodes present in the system: "))
    node_coords = []
    
    for i in range(node_num):        
        print(f"Enter coordinates for Node {i}:")
        
        x = float(input(f"  x-coordinate for Node {i}: "))
        y = float(input(f"  y-coordinate for Node {i}: "))
        z = float(input(f"  z-coordinate for Node {i}: "))

        print("Select node constraint type: \n 1 - Fixed Node \n 2 - Pinned Node \n 3 - Free Node")
        node_type = int(input(f" Enter the number corresponding to constraint type for Node {i}: "))

        if node_type == 1:
            constraints = [True, True, True, True, True, True]
        elif node_type == 2:
            constraints = [True, True, True, False, False, False]
        else:
            constraints = [False, False, False, False, False, False]
        
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

    for i in moment_nodes:
        if 0 <= i < len(node_pos):
            applied_forces[i][3:] = M[:]

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
    
    element_num = int(input("Enter the number of elements present in the system: "))
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
        
        element = Element(nodes[node_1_index], nodes[node_2_index], E, nu, A, I_y, I_z, J, local_z)
        element_connect.append(element)
        
    return element_connect

elements = define_elements(node_pos)

K_global = assemble_global_stiffness_matrix(node_pos, elements)

print("Defined Elements and Connected Nodes:")
for i, element in enumerate(elements, start=1):
    node_1, node_2 = element.nodes()
    length = element.length()
    print(f"Element {i}: connects Node({node_1.coordinates()}) and Node({node_2.coordinates()}) with length {length: .1f}")

print("Global Stiffness Matrix: ", K_global)


F_global = np.zeros(len(node_pos) * 6)

for node_idx, load in applied_loads.items():
    dof_index = node_idx * 6
    F_global[dof_index:dof_index+6] = load  # Assign applied forces and moments

K_mod, F_mod, constrained_DOFs = apply_boundary_conditions(K_global, F_global, node_pos)

# Solve for displacements
U_mod = solve_displacements(K_mod, F_mod)

# Reconstruct full displacement vector (including constrained DOFs)
U_full = np.zeros(len(node_pos) * 6)
free_dof_indices = [i for i in range(len(node_pos) * 6) if i not in constrained_DOFs]

for i, dof in enumerate(free_dof_indices):
    U_full[dof] = U_mod[i]

# Compute reaction forces
reaction_forces = compute_reaction_forces(K_global, U_full, constrained_DOFs)

# Display results
print("\nNodal Displacements:")
for i in range(len(node_pos)):
    print(f"Node {i}: Displacements = {U_full[i*6:i*6+3]}, Rotations = {U_full[i*6+3:i*6+6]}")

print("\nReaction Forces at Constrained DOFs:")
for dof in constrained_DOFs:
    print(f"DOF {dof}: {reaction_forces[dof]}")
