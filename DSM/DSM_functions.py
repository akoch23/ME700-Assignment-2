import math

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
    def __init__(self, node_1, node_2):
        
        self.node_1 = node_1
        self.node_2 = node_2
        
    def __repr__(self):
        return f"Element({self.node1}, {self.node2})"
    
    def nodes(self):
        return [self.node_1, self.node_2]
    
    def length(self):
        x_1, y_1, z_1 = self.node_1.coordinates()
        x_2, y_2, z_2 = self.node_2.coordinates()
        
        length = math.sqrt((x_2 - x_1)**2 + (y_2 - y_1)**2 + (z_2 - z_1)**2)
        return length

    
def node_coordinates():
    node_num = int(input("Enter the number of nodes present in the system: "))
    
    node_pos = []
    
    for i in range(node_num):        
        print(f"Enter coordinates for Node {i}:")
        
        x = float(input(f"  x-coordinate for Node {i}: "))
        y = float(input(f"  y-coordinate for Node {i}: "))
        z = float(input(f"  z-coordinate for Node {i}: "))
        
        new_node = Node(x, y, z)
        node_pos.append(new_node)
        
    return node_pos

node_pos = node_coordinates()

print("Defined Nodal Positions:")
for i, node in enumerate(node_pos, start=1):
    print(f"Node {i}: {node.coordinates()}")


def define_elements(nodes):
    element_num = int(input("Enter the number of elements present in the system: "))
    elements = []
    
    for i in range(element_num):
        print(f"Enter the nodes that serve as the end points for Element {i}:")
        
        node_1_index = int(input(f" Enter the index of Node 1 for Element {i}: "))
        node_2_index = int(input(f" Enter the index of Node 2 for Element {i}: "))
        
        element = Element(nodes[node_1_index], nodes[node_2_index])
        elements.append(element)
        
    return elements

elements = define_elements(node_pos)

print("Defined Elements and Connected Nodes:")
for i, element in enumerate(elements, start=1):
    node_1, node_2 = element.nodes()
    length = element.length()
    print(f"Element {i}: connects Node({node_1.coordinates()}) and Node({node_2.coordinates()}) with length {length: .1f}")



