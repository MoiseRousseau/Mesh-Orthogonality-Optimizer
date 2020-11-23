#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# Author : Moise Rousseau (2020), email at moise.rousseau@polymtl.ca

import numpy as np
import scipy.optimize as spopt
import multiprocessing
#import mpi4py as MPI #here store the element in each proc and make all collective op.

# in all the below
# idvertices is the id of the vertices as entered by the user
# ivertices is the id of the vertices in the structure in Ortho_Opt

from face_class import *


class Ortho_Opt:
  def __init__(self):
    self.penalizing_power = 1
    # input
    self.vertices = None
    self.elements = None
    self.elements_type = None
    self.elements_format = None
    self.initialized = False
    # optimization parameter
    self.numerical_derivative = False
    self.default_perturbation = 1e-6
    self.opt_algo_for_scipy = "BFGS"
    self.maxiter = 30
    self.stoptol = 1e-6
    # computed by Orth_Opt for work
    self.cell_center = None
    self.face_normal = None
    self.face_area = None
    self.face = None #contain the face instance
    self.fixed_vertices = set([0])
    
    #interesting array
    self.cost_function_value = None
    self.derivative = None
    self.i_iteration= 0
    return
    
  def set_penalizing_power(self, p):
    if p <= 0: 
      print("The penalizing power must be positive")
      print("Set to the default penalizing power (1)")
      self.penalizing_power = 1
      return
    self.penalizing_power = p
    return
  
  def use_numerical_derivative(self, x=True):
    self.numerical_derivative = x
    return
  
  def set_perturbation_numerical_derivative(self, x):
    self.default_perturbation = x
    return
    
  def set_optimisation_algorithm(self, algo):
    self.opt_algo_for_scipy = algo
    return
    
  def set_optimization_stop_tolerance(self,x):
    self.stoptol = x
    return
    
  def set_optimization_max_iteration(self,x):
    self.maxiter = x
    return
  
   
  ### OBJECTIVE FUNCTION ###
  def current_face_errors(self):
    if not self.initialized:
      self._initialize_all_()
      self.initialized = True
    self._compute_cell_center_all_()
    self._compute_cell_centers_vector_()
    self._compute_face_normal_and_area_()
    for iface, face in enumerate(self.face):
      error = face.compute_error()
      self.cost_function_value[iface] = face.weight * error ** self.penalizing_power
                                   
    bad_oriented_vector = self.cost_function_value > 1.
    if np.sum(bad_oriented_vector):
      for i,test in enumerate(bad_oriented_vector):
        if test: #bad oriented
          face = self.face[i]
          id_up = face.element_id_up
          face.element_id_up = face.element_id_dn
          face.element_id_dn = id_up
          v_up = face.vertice_up
          face.vertice_up = face.vertice_dn
          face.vertice_dn = v_up
          face.cell_centers_vector *= -1.
      #self._compute_cell_centers_vector_()
      for iface, face in enumerate(self.face):
        error = face.compute_error()
        self.cost_function_value[iface] = face.weight * error ** self.penalizing_power  
    return np.copy(self.cost_function_value)
  
  def current_cost_function(self):
    self.current_face_errors()
    return np.sum(self.cost_function_value)
  
  def current_face_orthogonality(self):
    self.current_face_errors()
    res = np.zeros(len(self.face), dtype='f8')
    for iface, face in enumerate(self.face):
      res[iface] = face.error
    return res
      
  
  def _compute_cell_center_all_(self):
    for ielem, elem in enumerate(self.elements):
      self.cell_center[ielem] = self._compute_cell_center_(ielem, elem)
    return
    
  def _compute_cell_center_(self, ielem, elem=None):
    if elem is None: elem = self.elements[ielem]
    elem = elem[elem != 0]
    nvertices = len(elem)
    return np.sum(self.vertices[elem-1],axis=0)/nvertices
  
  def _compute_cell_centers_vector_(self):
    #vector directed from self.face[-2] to self.face[-1]
    for face in self.face:
      id_up, id_dn = face.element_id_up, face.element_id_dn
      face.cell_centers_vector = self.cell_center[id_up-1] - self.cell_center[id_dn-1]
      face.cell_centers_vector_norm = np.linalg.norm(face.cell_centers_vector)
      face.cell_centers_vector /= face.cell_centers_vector_norm 
    return
    
  def _compute_face_normal_and_area_(self):
    for face in self.face:
      u = self.vertices[face.vertices[1]-1] - self.vertices[face.vertices[0]-1]
      v = self.vertices[face.vertices[2]-1] - self.vertices[face.vertices[1]-1]
      face.normal = np.cross(u,v)
      face.area = np.linalg.norm(face.normal)
      face.normal /= face.area
      if face.type == 3: face.area /= 2
    return
    
  
    
  ### DERIVATIVE ###
  def current_derivative(self):
    if not self.initialized:
      self._initialize_all_()
      self.initialized = True
    self.current_face_errors()
    self.derivative[:] = 0.
    if self.numerical_derivative:
      for iface,face in enumerate(self.face):
        for ivertice in face.vertices: #face vertices
          if ivertice in self.fixed_vertices: continue 
          self._derivative_finite_difference_(ivertice, iface, face)
        if face.vertice_up not in self.fixed_vertices:
          self._derivative_finite_difference_(face.vertice_up, iface, face)
        if face.vertice_dn not in self.fixed_vertices:
          self._derivative_finite_difference_(face.vertice_dn, iface, face)
    else:
      for iface,face in enumerate(self.face):
        element_id1_type = self.elements_type[face.element_id_dn-1]
        element_id2_type = self.elements_type[face.element_id_up-1]
        if element_id1_type == 4 and element_id2_type == 4: #case 1: tet / tet
          self._derivative_case_1_(face)
        elif element_id1_type == 4 and element_id2_type == 5: #case 2: tet / pyr
          pass
        elif (element_id1_type == 5 and element_id2_type == 5
              and face[3] == 0): #case 3: tri pyr / pyr
          pass
        elif (element_id1_type == 5 and element_id2_type == 5
              and face[3] != 0): #case 4: quad pyr / pyr
          pass
        else:
          #TODO what to do with prisms and hex
          pass
      #TODO see what happen...
      self.derivative *= -1.
    return self.derivative
  
    
  def _derivative_finite_difference_(self, idvertice, iface, face=None):
    if face is None: face = self.face[iface]
    for idir in [0,1,2]: #X,Y,Z
      current_error = self.cost_function_value[iface]
      #perturbate the position
      self.vertices[idvertice-1,idir] += self.default_perturbation
      #compute cell centers vector
      E = self._compute_cell_center_(face.element_id_dn-1)
      F = self._compute_cell_center_(face.element_id_up-1)
      r_f = F-E
      r_f /= np.linalg.norm(r_f)
      #normal
      u = self.vertices[face.vertices[1]-1] - self.vertices[face.vertices[0]-1]
      v = self.vertices[face.vertices[2]-1] - self.vertices[face.vertices[1]-1]
      n_f = np.cross(u,v)
      n_f /= np.linalg.norm(n_f)
      #new error
      pertur_error = (1 - r_f.dot(n_f))
      pertur_error = face.weight * pertur_error ** self.penalizing_power
      #compute derivative
      deriv = (pertur_error - current_error) / self.default_perturbation
      self.derivative[idvertice-1,idir] += deriv
      #restore position
      self.vertices[idvertice-1,idir] -= self.default_perturbation
    return
    
  def _derivative_case_1_(self, face):
    #treat first vertice in A position
    for i, ivertice in enumerate(face.vertices):
      if ivertice in self.fixed_vertices: continue
      C = face.vertices[i-1]
      if i == 2: B = face.vertices[0]
      else: B = face.vertices[i+1]
      vec = (face.cell_centers_vector - (1-face.error) * face.normal)
      deriv = np.cross(vec, self.vertices[C-1] - self.vertices[B-1])
      self.derivative[ivertice-1] += deriv / (2. * face.area) 
    #F position
    vec = (face.normal - (1-face.error) * face.cell_centers_vector)
    deriv = 0.25 * vec / face.cell_centers_vector_norm * \
            face.weight * face.error ** (self.penalizing_power-1)
    if face.vertice_dn not in self.fixed_vertices:
      self.derivative[face.vertice_dn-1] -= deriv #TODO check sign
    #E position
    if face.vertice_up not in self.fixed_vertices:
      self.derivative[face.vertice_up-1] += deriv
    return
    
  
  
  ### OPTIMIZATION ###
  def _wrapper_cost_funtion_scipy_minimize_(self, x):
    self.i_iteration += 1
    self.load_vertices(np.reshape(x, (int(len(x)/3),3)))
    res = self.current_cost_function()
    print(f"Iteration {self.i_iteration} / Current cost function: {res}")
    return res
  
  def _wrapper_derivative_scipy_minimize_(self, x):
    #self.load_vertices(np.reshape(x, (int(len(x)/3),3)))
    res = self.current_derivative()
    #print(res[450:])
    return np.ravel(res[:])
    
    
  def optimize(self, algo=None, maxiter=None, tol=None):
    if not algo: algo = self.opt_algo_for_scipy
    if not maxiter: maxiter = self.maxiter
    if not tol: tol=self.stoptol
    x0 = np.copy(np.ravel(self.vertices))
    try:
      res = spopt.minimize(self._wrapper_cost_funtion_scipy_minimize_, x0, 
                         jac=self._wrapper_derivative_scipy_minimize_,
                         tol = tol,
                         options={"maxiter":maxiter})
    except KeyboardInterrupt:
      return
    print("\nOptimization terminated. Output from SciPy:")
    print(res)
    self.load_vertices(np.reshape(res.x, (int(len(res.x)/3),3)))
    return
  
  def optimize_nlopt(self, tol=None, maxiter=None):
    if not tol: tol=self.stoptol
    if not maxiter: maxiter = self.maxiter
    optimizer = nlopt.opt(nlopt.LD_LBFGS, len(np.ravel(self.vertices)))
    optimizer.set_min_objective(self._wrapper_nlopt_)
    optimizer.set_ftol_rel(tol)
    optimizer.set_maxeval(maxiter)
    x0 = np.ravel(self.vertices)
    res = optimizer.optimize(x0)
    self.load_vertices(np.reshape(res, (int(len(res)/3),3)))
    return
  
  
  ### LOAD VERTICES ###
  def load_vertices(self, nparray):
    # A Numpy array of size (n,3)
    self.vertices = nparray
    return
    
  def define_fixed_vertices(self, fixed):
    self.fixed_vertices = np.array(fixed,dtype='i8')
    return
  
  def load_vertices_from_file(self, f, f_format=None):
    if f_format is None:
      print("Deduce vertices file format from extension")
      f_format = f.split('.')[-1]
    f_format = f_format.lower()
    self.vertices_format = f_format
    if f_format == "node": #tetgen tetra
      self.vertices = np.genfromtxt(f,usecols = (1,2,3), comments="#")
    f_format = f_format.lower()
    return
    
  ### SAVE VERTICES ###
  def write_vertices(self, fout, f_format=None):
    if f_format is None:
      print("Deduce vertices file format from extension")
      f_format = fout.split('.')[-1]
    if f_format == "numpy":
      np.savetxt(fout, self.vertices)
    return
      
  
  
  ### LOAD ELEMENTS ###
  def load_elements(self, nparray):
    #np array of max size 8 and containing 0 for ignoring vertices
    #Be carefull in element ordering
    self.elements = nparray
    return 
  
  def load_elements_from_file(self, f, f_format=None):
    if f_format is None:
      print("Deduce elements file format from extension")
      f_format = f.split('.')[-1]
    f_format = f_format.lower()
    self.elements_format = f_format
    if f_format == "ele": #tetgen tetra
      self.elements = np.genfromtxt(f,usecols = (1,2,3,4), comments="#", dtype='i8')
      self.elements_type = np.ones(len(self.elements), dtype="i2") * 4 #tet
      pass
    elif f_format == "med": #salome
      pass
    #TODO order here the numbering ??
    return
    
  def _create_face_info_(self):
    #from the input, we need to identify all the face in the mesh
    # therefore, for each element, store the faces and the two element id
    #
    # for each face, store the 3/4 considered face vertices and the two
    # element id which share the face
    # and also the 2 others point (E and F)
    print("\nBuilding internal faces")
    nfaces = 0 #number of face
    temp_dict = {} #hold element shared and other vertices (E, F)
    tet_faces = [np.array([0,1,2],dtype='i2'), np.array([1,2,3],dtype='i2'),
                 np.array([0,2,3],dtype='i2'), np.array([0,1,3],dtype='i2')]
    pyr_faces = [np.array([0,1,2],dtype='i2'), np.array([1,2,3],dtype='i2'),
                 np.array([0,2,3],dtype='i2'), np.array([0,1,3],dtype='i2')]

    for i,elem in enumerate(self.elements):
      elem_type = np.sum(elem != 0) #count number of vertices in the element
      if elem_type == 4: #Tet
        for ifacetest, f in enumerate(tet_faces):
          face_vertices = elem[f].tolist()
          min_index = face_vertices.index(min(face_vertices))
          if (min_index != len(face_vertices) -1 and 
              face_vertices[min_index+1] < face_vertices[min_index-1]):
            face_vertices = tuple(face_vertices[min_index:] + face_vertices[:min_index])
          elif (min_index == len(face_vertices) - 1 and 
               face_vertices[0] < face_vertices[min_index-1]):
            face_vertices = tuple(face_vertices[min_index:] + face_vertices[:min_index])
          else:
            face_vertices.reverse()
            min_index = len(face_vertices) - min_index -1
            face_vertices = tuple(face_vertices[min_index:] + face_vertices[:min_index])
          #add it to the temp dict
          try: 
            temp_dict[face_vertices][0:2] = [i+1, elem[ifacetest-1]] #id_up, vertice_up
            nfaces += 1
          except: 
            temp_dict[face_vertices] = [0, 0, i+1, elem[ifacetest-1]]
      elif elem_type == 5: #pyramid
        pass
      elif elem_type == 6: #prisms
        pass
      elif elem_type == 8: #Hex
        pass
    
    #now we have all the face and have the id of both element sharing the face
    #self.face is a list of the face instance
    print(f"{nfaces} internal faces detected")
    self.face = [Face() for x in range(nfaces)]
    count = 0
    for face_vertices,ids in temp_dict.items():
      if not ids[0]: #face is boundary face (not shared)
        for x in face_vertices: self.fixed_vertices.add(x) #boudnary point (1+ based id)
      elif len(face_vertices) == 3:
        face = self.face[count]
        face.populate_face(ids[0], ids[2], face_vertices, ids[1], ids[3])
        count += 1
      elif len(face_vertices) == 4:
        for x in face_vertices: self.fixed_vertices.add(x)
        self.face[count,:4] = face_vertices
        self.face[count,4:] = ids
        count += 1
    print(f"{len(self.fixed_vertices)-1} boundary points in the mesh")
    return
  
  
  
  ### INITIALIZATION FUNCTION ###
  def _initialize_all_(self):
    if not self.face: self._create_face_info_() #can also be user supplied
    self.face_normal = np.zeros((len(self.face),3),dtype='f8')
    self.face_area = np.zeros(len(self.face),dtype='f8')
    self.cell_center = np.zeros((len(self.elements),3), dtype='f8')
    self.cost_function_value = np.zeros(len(self.face), dtype='f8')
    self.derivative = np.zeros((len(self.vertices),3),dtype='f8')
    return
    
