# Step 1: Define Nodes
syntax: 
nodes = {node ID: np.array([x, y, z])

Example:
nodes = {  
    0: np.array([0, 0, 0]),  
    1: np.array([-5, 1, 10]),  
    2: np.array([-1, 5, 13]),  
    3: np.array([-3,7,11]),  
    4: np.array([6,9,5])  
}

# Step 2: Define Element Properties
syntax: 
section_props_element_# = {"E": #, "nu": #, "A": #, "Iz": #, "Iy": #, "J": #, "local_z": np.array([#, #, #])}

Example:
section_props_element_1 = {
    "E": 500, "nu": 0.3, "A": np.pi, "Iz": np.pi/4, "Iy": np.pi/4, "J": np.pi/2, "local_z": np.array([0.0, 0.0, 1.0])
}

section_props_element_2 = {
    "E": 1000, "nu": 0.3, "A": 0.5, "Iz": 0.04167, "Iy": 0.01042, "J": 0.02861, "local_z": np.array([1.0, 0.0, 0.0])
}

section_props_element_3 = {
    "E": 210e9, "nu": 0.3, "A": 0.01, "Iz": 8.33e-6, "Iy": 8.33e-6, "J": 1.67e-5, "local_z": np.array([0.0, 0.0, 1.0])
}

Use less/add more elements as needed

# Step 3: Define Elements
syntax: 
elements = [
(nodeID1, nodeID2, section_props_element_1),  # first element  
(nodeID1, nodeID2, section_props_element_2)   # second element
]

Example:
elements = [
    (0, 1, section_props_element_1),
    (1, 2, section_props_element_1),
    (2, 3, section_props_element_1),
    (2, 4, section_props_element_1)
]

# Step 4: Apply Nodal Loads
syntax:
loads = {
  nodeID: np.array([Fx, Fy, Fz, Mx, My, Mz])
}

Example:
loads = {
    1: np.array([0.1, -0.05, -0.075, 0, 0, 0]),
    2: np.array([0, 0, 0, 0.5, -0.1, 0.3])
}

# Step 5: Apply Boundary Conditions (Constrain Nodal DOFs)
syntax:
supports = {
  nodeID: [True, True, True, True, True, True], # Completely fixed
  nodeID: [False, True, True, False, False, True] # Motion allowed in x direction and rotation allowed along the x and y axis
}

Example:
supports = {
    0: [False, False, True, False, False, False],
    4: [True, True, True, False, False, False],
    3: [True, True, True, True, True, True]
}

# Step 6: Use Direct Stiffness Method Solver to solve for Nodal Displacements/Rotations and Reaction Forces/Moments

Initialize solver:
solver = dsm.Frame3DSolver(nodes, elements, loads, supports)
displacements, reactions = solver.solve()

disp_matrix = displacements.reshape((-1, 6))
reac_matrix = reactions.reshape((-1, 6))

Create a dictionary for displacements and reactions
disp_dict = {node: disp_matrix[i] for i, node in enumerate(nodes)}
react_dict = {node: reac_matrix[i] for i, node in enumerate(nodes)}

To output results:
print("Nodal Displacements and Rotations:")
for node, disp in disp_dict.items():
  print(f"Node {node}: [u: {disp[0]:.10f}, v: {disp[1]:.10f}, w: {disp[2]:.10f}, "
        f"rot_x: {disp[3]:.10f}, rot_y: {disp[4]:.10f}, rot_z: {disp[5]:.10f}]")
    
print("\nReaction Forces and Moments at Supports:")
for node, react in react_dict.items():
  Only display reactions for nodes with boundary conditions
  print(f"Node {node}: [Fx: {react[0]:.10f}, Fy: {react[1]:.10f}, Fz: {react[2]:.10f}, "
        f"Mx: {react[3]:.10f}, My: {react[4]:.10f}, Mz: {react[5]:.10f}]")

Example Output:
Nodal Displacements and Rotations:
Node 0: [u: 0.1629742239, v: 0.0675373057, w: 0.0000000000, rot_x: 0.0038603044, rot_y: -0.0097767200, rot_z: 0.0099770437]
Node 1: [u: 0.0568382275, v: -0.0212726063, w: -0.0442346868, rot_x: 0.0039556640, rot_y: -0.0092999220, rot_z: 0.0099770437]
Node 2: [u: 0.0010444257, v: 0.0010905207, w: 0.0003463178, rot_x: 0.0031355569, rot_y: -0.0040054958, rot_z: 0.0051427325]
Node 3: [u: 0.0000000000, v: 0.0000000000, w: 0.0000000000, rot_x: 0.0000000000, rot_y: 0.0000000000, rot_z: 0.0000000000]
Node 4: [u: 0.0000000000, v: 0.0000000000, w: 0.0000000000, rot_x: -0.0045516627, rot_y: 0.0004901880, rot_z: 0.0006642572]

Reaction Forces and Moments at Supports:
Node 0: [Fx: 0.0000000000, Fy: 0.0000000000, Fz: 0.0066721997, Mx: 0.0000000000, My: 0.0000000000, Mz: -0.0000000000]
Node 1: [Fx: -0.0000000000, Fy: -0.0000000000, Fz: 0.0000000000, Mx: 0.0000000000, My: 0.0000000000, Mz: -0.0000000000]
Node 2: [Fx: 0.0000000000, Fy: 0.0000000000, Fz: -0.0000000000, Mx: -0.0000000000, My: -0.0000000000, Mz: 0.0000000000]
Node 3: [Fx: -0.0235127129, Fy: 0.1379482485, Fz: 0.0253249828, Mx: -0.4116107460, My: 0.2981182342, Mz: -0.3614403375]
Node 4: [Fx: -0.0764872871, Fy: -0.0879482485, Fz: 0.0430028175, Mx: 0.0000000000, My: -0.0000000000, Mz: 0.0000000000]

# Step 7: Plotting Reaction Forces/Momenets for Individual Elements and Deformed Structure

Initialize plot function (Elements):
internal_forces = solver.compute_internal_forces_and_moments(displacements)
solver.plot_internal_forces_and_moments(internal_forces)

Initialize plot function (Structure)
solver.plot_deformed_shape(displacements, scale=25)

# Step 8: Computing Critical Load Factor and Plooting Buckling Mode Shapes

Define nodes and elements for buckling analysis
nodes_ecls = {
    0: np.array([0.0, 0.0, 0.0]),
    1: np.array([30.0, 40.0, 0.0])
}
User can input element properties as demonstrated below. If I_rho is not provided then function will take I_rho = J as default.

elements_ecls = [
    (0, 1, {"E": 1000, "nu": 0.3, "A": np.pi, "Iz": np.pi/4, "Iy": np.pi/4, "J": np.pi/2, "I_rho": np.pi/2})
]

loads_ecls = {
    1: np.array([-3/5, -4/5, 0, 0, 0, 0])
}

supports_ecls = {
    0: [True, True, True, True, True, True]
}

frame_solver_ecla = dsm.Frame3DSolver(nodes_ecls, elements_ecls, loads_ecls, supports_ecls)

Solve for buckling modes WITH and WITHOUT the interaction terms
for use_interaction in [True, False]:
    solver_type = "Without Interaction Terms" if not use_interaction else "With Interaction Terms"
    print(f"Solving for {solver_type}")

    ecl_solver = ecls.ElasticCriticalLoadSolver(frame_solver_ecla, use_interaction_terms=use_interaction)
    eigenvalues, eigenvectors = ecl_solver.solve_eigenvalue_problem()

    mode_shape = eigenvectors[:, 0]

    # Plot the first buckling mode using Hermite shape functions
    print("Critical Load Factors:", eigenvalues)
    print(f"Lowest Critical Load Factor: {np.min(eigenvalues)}")
    ecls.plot_buckling_mode(ecl_solver.frame_solver, mode_shape, scale_factor=5)
