#!/usr/bin/env python
# coding: utf-8

# In[ ]:


node(self, u, v, w, theta_x, theta_y, theta_z):
    self.u = u
    self.v = v
    self.w = w
    self.theta_x = theta_x
    self.theta_y = theta_y
    self.theta_z = theta_z
    
    def node_coords(self):
        return [self.u, self.v, self.w, self.theta_x, self.theta_y, self.theta_z]
    
element(self, node_1, node_2):
    self.node_1 = node_1
    self.node_2 = node_2
    
    def element_coords(self):
        return [self.node_1.node_coords(), self.node_2.node_coords()]
    
def stiffness_matrix ():

