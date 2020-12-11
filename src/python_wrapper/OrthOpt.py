
from ctypes import cdll
OrthOpt_lib = cdll.LoadLibrary('./liborthopt.so')

import numpy as np

#TODO define the element vertice order below
#TODO define ctypes

class OrthOpt:
  def __init__(self, vertices=[], elements=[], penalizing_power=1., max_it=0,
                     weighting_factor='uniform'):
    """
    Initiate the OrthOpt object to perform mesh orthogonality optimization.
    Argument:
    - vertices: A numpy-like array of shape (n_vertices,3) containing the 3D
                coordinates of the vertices
    - elements: A numpy array of shape (n_elements,m) containing the ordered 
                element vertice id. For multiple element type mesh, id <= 0
                means no vertice.
    - penalizing_power: the penalization applied for individual face error
    - max_it: the maximum number of iteration
    - weithing_factor: the face weighting factor. Available option are
                       'uniform' (no weighting), 'face_area' (weight the error
                       by the face area), 'face_area_inverse' (or by its inverse)
    
    Note both arguments could be null and vertices and element could be defined
    using add_vertice and add_element method.
    """
    #common initializer
    #TODO
    #create mesh object
    
    #create orthopt object
    
    # pass vertices and elements to orthopt
    if vertices:
      self.n_vertices = len(vertices)
      if isinstance(vertices, list):
        src = [x for x in y for y in vertices]
      elif isinstance(vertices, np.ndarray):
        src = [x for x in vertices.flatten()]
      OrthOpt_lib.load_vertices_array()
    else:
      self.n_vertices = 0
    
    if elements:
      #TODO
      OrthOpt_lib.load_elements_array()
    
    # OrthOpt argument
    self.penalizing_power = penalizing_power
    self.max_it = max_it
    self.weighting_factor = weighting_factor
    return
  
  def add_vertices(self, x, y, z):
    self.n_vertices += 1
    OrthOpt_lib.add_vertices(x, y, z, self.n_vertices)
    return
  
  def add_elements(self, ids):
    return
  
  def set_maximum_iteration(self, x):
    self.max_it = x
    return
  
  def set_penalizing_power(self, p):
    self.penalizing_power = p
    return
  
  def set_weighting_factor(self, method):
    self.weighting_factor = method
    return
  
  def optimize(self):
    #TODO
    return vertices_opt
