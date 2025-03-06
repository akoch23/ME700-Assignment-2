import numpy as np
# Given Node positions, 3D Force, 3D Moment, b, h, E, v, J, E0, E1, A, I_y, I_z, I_rho

def node_coordinates():
  # Asking user to input number of nodes in the system
  node_num = int(input("Enter the numbe of nodes in the system: ")

  # List for storing coordinate sets
  node_coords = []
                 
  for i in range(node_num):
    print(f"Enter the coordinates for Node {i+1}:")
    x = float(input("x coordinate: ")
    y = float(input("y coordinate: ")
    z = float(input("z coordinate: ")

    node_coords.append((x,y,z))

  print("\nNode coordinates:")
  for i, coord in enumerate(node_coords, start=1):
        print(f"Node {i}: x = {coord[0]}, y = {coord[1]}, z = {coord[2]}")
        
  return node_coords

node_coordinates
  

