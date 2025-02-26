
node(self, u, v, w, theta_x, theta_y, theta_z, f_x, f_y, f_z, mu_x, mu_y, mu_z)):
    self.u = u
    self.v = v
    self.w = w
    self.theta_x = theta_x
    self.theta_y = theta_y
    self.theta_z = theta_z
    self.f_x = f_x
    self.f_y = f_y
    self.f_z = f_z
    self.mu_x = mu_x
    self.mu_y = mu_y
    self.mu_z = mu_z
    
    def node_dof(self):
        return [self.u, self.v, self.w, self.theta_x, self.theta_y, self.theta_z]

    def node_loads(self):
        return [self.f_x, self.f_y, self.f_z, self.mu_x, self.mu_y, self.mu_z]
    
element(self, node_1, node_2):
    self.node_1 = node_1
    self.node_2 = node_2
    
    def displacement_rot_vector(self):
        return [self.node_1.node_dof(), self.node_2.node_dof()]
        
    def load_vector(self):
        return [self.node_1.node_loads(), self.node_2.node_loads()]

def stiffness_matrix(F, k, E, A, L, u):
    
    
    
