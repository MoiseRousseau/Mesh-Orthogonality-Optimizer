
import ctypes
OrthOpt_lib = ctypes.cdll.LoadLibrary('./liborthopt.so')

import numpy as np

#TODO define the element vertice order below

class OrthOpt:
  def __init__(self, vertices=None, elements=None, penalizing_power=1., max_it=0,
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
    #parse input
    self.penalizing_power = penalizing_power
    self.max_it = max_it
    self.weighting_factor = weighting_factor
    
    #create mesh object
    self.mesh = OrthOpt_lib.Mesh_init()
    if vertices is not None:
      self.n_vertices = len(vertices)
      if isinstance(vertices, list):
        _vertices = np.array([x for x in y for y in vertices], dtype='f8')
      elif isinstance(vertices, np.ndarray):
        _vertices = vertices.flatten()
      print("here")
      OrthOpt_lib.Mesh_load_vertices_array(self.mesh, _vertices.ctypes.data, 
                                           ctypes.c_int(self.n_vertices))
      print("here")
    else:
      self.n_vertices = 0
    if elements is not None:
      self.n_elements = len(elements)
      if isinstance(elements, list):
        _elements = np.array([x for x in y for y in elements], dtype='i8')
        _types = np.array([len(x) for x in elements], dtype='i8')
      elif isinstance(elements, np.ndarray):
        _elements = elements.flatten()
        _types = np.ones(self.n_elements, dtype='i8')*len(_elements)/self.n_elements
      OrthOpt_lib.Mesh_load_elements_array(self.mesh, _elements.ctypes.data,
                                           _types.ctypes.data, self.n_elements)
    else:
      self.n_elements = 0
    
    #create orthopt object
    
    
    return
  
  #MESH METHOD
  def add_vertices(self, x, y, z):
    self.n_vertices += 1
    OrthOpt_lib.Mesh_add_vertices(self.mesh, ctypes.c_double(x), 
                             ctypes.c_double(y), ctypes.c_double(z), 
                             self.n_vertices)
    return
  def add_elements(self, ids):
    return
  
  #OPTIMIZER PARAMETER
  def set_maximum_iteration(self, x):
    self.max_it = x
    return
  def set_penalizing_power(self, p):
    self.penalizing_power = p
    return
  def set_weighting_factor(self, method):
    self.weighting_factor = method
    return
  
  #OPTIMIZER
  def get_face_error(self):
    return
  def get_derivative(self):
    return
  def optimize(self):
    #TODO
    return vertices_opt
