import math
import numpy as np

# Given Node positions, 3D Force, 3D Moment, b, h, E, v, J, E0, E1, A, I_y, I_z, I_rho

class Node:
    def __init__(self, x, y, z,):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"Node({self.x}, {self.y}, {self.z})"    

    def coordinates(self):
        return [self.x, self.y, self.z]
    
    
class Element:
    def __init__(self, node_1, node_2, E, nu, A, I_y, I_z, I_rho, J):
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


def node_coordinates():
    node_num = int(input("Enter the number of nodes present in the system: "))
    node_coords = []
    
    for i in range(node_num):        
        print(f"Enter coordinates for Node {i}:")
        
        x = float(input(f"  x-coordinate for Node {i}: "))
        y = float(input(f"  y-coordinate for Node {i}: "))
        z = float(input(f"  z-coordinate for Node {i}: "))
        
        new_node = Node(x, y, z)
        node_coords.append(new_node)
        
    return node_coords

node_pos = node_coordinates()

print("Defined Nodal Positions:")
for i, node in enumerate(node_pos, start=1):
    print(f"Node {i}: {node.coordinates()}")


def define_elements(nodes):
    element_num = int(input("Enter the number of elements present in the system: "))
    element_connect = []
    
    for i in range(element_num):
        print(f"Enter the nodes that serve as the end points for Element {i}:")
        
        node_1_index = int(input(f" Enter the index of Node 1 for Element {i}: "))
        node_2_index = int(input(f" Enter the index of Node 2 for Element {i}: "))
        
        element = Element(nodes[node_1_index], nodes[node_2_index])
        element_connect.append(element)
        
    return element_connect

elements = define_elements(node_pos)

print("Defined Elements and Connected Nodes:")
for i, element in enumerate(elements, start=1):
    node_1, node_2 = element.nodes()
    length = element.length()
    print(f"Element {i}: connects Node({node_1.coordinates()}) and Node({node_2.coordinates()}) with length {length: .1f}")



