import math
import numpy as np

# Given Node positions, 3D Force, 3D Moment, b, h, E, v, J, E0, E1, A, I_y, I_z, I_rho

class Node:
    def __init__(self, x, y, z, constraints=None, forces=None):
        self.x = x
        self.y = y
        self.z = z
        self.constraints = constraints if constraints else [False] * 6
        self.force = force if force else [0] * 6

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
        self.I_z I_z
        self.I_rho = I_rho
        self.J = J
        
    def __repr__(self):
        return f"Element({self.node1}, {self.node2})"
    
    def nodes(self):
        return [self.node_1, self.node_2]
    
    def length(self):
        x_1, y_1, z_1 = self.node_1.coordinates()
        x_2, y_2, z_2 = self.node_2.coordinates()
        
        length = math.sqrt((x_2 - x_1)**2 + (y_2 - y_1)**2 + (z_2 - z_1)**2)
        return length

    def stiffness_matrix(self):
        L = self.length()
        K_local = local_elastic_stiffness_matrix_3D_beam(self.E, self.A, self.I_y, self.I_z, self.J, L)

        rotation = rotation_matrix_3D(self.node_1.coordinates(), self.node_2.coordinates())
        transformation = transformation_matrix_3D(rotation)

        K_global = transformation.T @ k_local @ transformation
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
                K_global[dof_indices[i], dof_indices[j]], += K_e[i,j]
    return K_global
'''
def apply_boundary_conditions(K, F, nodes):
    constrained_DOFs = []
    for i, node in enumerate(nodes):
        for dof in range(6):
            if node.constraints[dof]:
                constrained_dofs.append(i * 6 + dof)

    K_mod = np.delete(K, constrained_DOFs, axis = 0)
    K_mod = np.delete(K_mod, constrained_DOFs, axis=1)
    F_mod = np.delete(F, constrained_DOFs, axis=0)

    return K_mod, F_mod, constrained_DOFs

def solve_displacements(K_mod, F_mod):
    return np.linalg.solve(K_mod, F_mod)

def compute_reaction_forces(K, U_full, constrained_DOFs):
    return K @ U_full
'''

def node_coordinates():
    node_num = int(input("Enter the number of nodes present in the system: "))
    node_coords = []
    
    for i in range(node_num):        
        print(f"Enter coordinates for Node {i}:")
        
        x = float(input(f"  x-coordinate for Node {i}: "))
        y = float(input(f"  y-coordinate for Node {i}: "))
        z = float(input(f"  z-coordinate for Node {i}: "))

        print("Select node constraint type: \n 1 - Fixed Node \n 2 - Pinned Node \n 3 - Free Node")
        node_type = int(input(" Enter the number corresponding to contraint type for Node {i}: "))

        if node_type == 1:
            constraints = [True, True, True, True, True, True]
        elif node_type == 2:
            constraints = [True, True, True, False, False, False]
        else:
            constraints = [False, False, False, False, False, False]
        
        new_node = Node(x, y, z)
        node_coords.append(new_node)
    
    return node_coords

def define_forces_and_moments():
    F_x = float(input(f" Normal force acting in x-direction: ") 
    F_y = float(input(f" Normal force acting in y-direction: ") 
    F_z = float(input(f" Normal force acting in z-direction: ") 

    F = [F_x, F_y, F_z]

    M_x = float(input(f" Moment force acting in x-direction: ") 
    M_y = float(input(f" Moment force acting in y-direction: ") 
    M_z = float(input(f" Moment force acting in z-direction: ") 

    M = [M_x, M_y, M_z]

    return F, M

def apply_forces_and_moments(nodes, F, M):
    
    force_input = input("Enter the node where Normal forces (F_x, F_y, F_z) are applied: ")
    if force_input.strip():
        node_idx = int(force_input)
        if node_idx < len(nodes):
            print(f"Force at Node {node_idx}: Fx = {F[node_idx][0]}, Fy = {F[node_idx][1]}, Fz = {F[node_idx][2]}")
            
    moment_input = input("Enter the node where Moment forces (M_x, M_y, M_z) are applied: ")
    if moment_input.strip():
        node_idx = int(moment_input)
        if node_idx < len(nodes):
            print(f"Moment at Node {node_idx}: Mx = {M[node_idx][0]}, My = {M[node_idx][1]}, Mz = {M[node_idx][2]}")

    return F, M

node_pos = node_coordinates()

F = np.zeros((len(nodes), 3))
M = np.zeros((len(nodes), 3))

F, M = define_forces_and_moments(nodes, F, M)


print("Defined Nodal Positions:")
for i, node in enumerate(node_pos, start=1):
    print(f"Node {i}: {node.coordinates()}")

print("\nApplied Forces (F):")
for i, force in enumerate(F):
    print(f"Node {i}: Fx = {force[0]}, Fy = {force[1]}, Fz = {force[2]}")

print("\nApplied Moments (M):")
for i, moment in enumerate(M):
    print(f"Node {i}: Mx = {moment[0]}, My = {moment[1]}, Mz = {moment[2]}")

def define_elements(nodes):
    element_num = int(input("Enter the number of elements present in the system: "))
    element_connect = []
    
    for i in range(element_num):
        print(f"Enter the nodes that serve as the end points for Element {i}:")
        
        node_1_index = int(input(f" Enter the index of Node 1 for Element {i}: "))
        node_2_index = int(input(f" Enter the index of Node 2 for Element {i}: "))

        E = Element.E
        nu = Element.nu
        A = Element.A
        I_y = Element.I_y
        I_z = Element.I_z
        J = Element.J
        
        element = Element(nodes[node_1_index], nodes[node_2_index], E, nu, A, I_y, I_z, J)
        element_connect.append(element)
        
    return element_connect

elements = define_elements(node_pos)

print("Defined Elements and Connected Nodes:")
for i, element in enumerate(elements, start=1):
    node_1, node_2 = element.nodes()
    length = element.length()
    print(f"Element {i}: connects Node({node_1.coordinates()}) and Node({node_2.coordinates()}) with length {length: .1f}")
print("Global Stiffness Matrix: "K_global)



