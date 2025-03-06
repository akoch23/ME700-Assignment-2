import numpy as np
# Given Node positions, 3D Force, 3D Moment, b, h, E, v, J, E0, E1, A, I_y, I_z, I_rho

def node_coordinates():
    node_num = int(input("Enter the numbe of nodes present in the system: "))
    
    node_pos = []
    
    for i in range(node_num):        
        print(f"Enter coordinates for Node {i}:")
        
        x = float(input(f"  x-coordinate for Node {i}: "))
        y = float(input(f"  y-coordinate for Node {i}: "))
        z = float(input(f"  z-coordinate for Node {i}: "))

        node_pos.append([x,y,z])
        
    return node_pos

node_pos = node_coordinates()

print("Defined Nodal Positions:")
for i, node in enumerate(node_pos, start=1):
    print(f"Node {i}: {node}")
  

  

