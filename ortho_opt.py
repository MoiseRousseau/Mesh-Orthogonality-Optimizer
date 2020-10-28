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
#import mpi4py as MPI #here store the element in each proc and make all collective op.

TWO_TET = 0 #triangle shared by two tetrahedra
TRI_TET_PYR = 1 #triangle shared by one tet and one pyramid
TRI_TWO_PYR = 2 #triangle shared by two pyramids
QUAD_TWO_PYR = 3 #quad shared by two pyramids


class Ortho_Opt:
  def __init__(self):
    self.penalizing_power = 1
    # input
    self.vertices = None
    self.elements = None
    self.elements_format = None
    self.initialized = False
    # computed by Orth_Opt
    self.cell_center = None
    self.cell_centers_vector = None
    self.face_normal_area = None #product of face normal and area
    self.face = None #store the 3/4 considered face vertices and the two element ids
    
    return
    
  def set_penalizing_power(self, p):
    if p <= 0: 
      print("The penalizing power must be positive")
      print("Set to the default penalizing power (1)")
      self.penalizing_power = 1
      return
    self.penalizing_power = p
    return
  
   
  ### OBJECTIVE FUNCTION ###
  def current_cost_function(self):
    if not self.initialized:
      self._initialize_all_()
      self.initialized = True
    self._compute_cell_center_()
    self._compute_cell_centers_vector_()
    self._compute_face_normal_area_()
    current = (self.cell_centers_vector[:,0] * self.face_normal_area[:,0] +
               self.cell_centers_vector[:,1] * self.face_normal_area[:,1] +
               self.cell_centers_vector[:,2] * self.face_normal_area[:,2]) ** self.penalizing_power
    return np.sum(current)
    
  def _compute_cell_center_(self):
    for i,elem in enumerate(self.elements):
      count = 0
      for v in elem:
        if not v: continue
        self.cell_center[i] += self.vertices[v-1,:]
        count += 1
      self.cell_center[i] /= count
    return
  
  def _compute_cell_centers_vector_(self):
    #vector directed from self.face[-2] to self.face[-1]
    self.cell_centers_vector[self.face[:,-2] != 0] = (self.cell_center[self.face[:,-1]-1,:] -
            self.cell_center[self.face[:,-2]-1,:]) [self.face[:,-2] != 0]
    #cell center vector on the boundary is 0
    return
    
  def _compute_face_normal_area_(self):
    u = self.vertices[self.face[:,1]-1] - self.vertices[self.face[:,0]-1]
    v = self.vertices[self.face[:,2]-1] - self.vertices[self.face[:,0]-1]
    self.face_normal_area[:] = np.cross(u,v,axis=1)
    self.face_normal_area[self.face[:,3] == 0] /= 2
    return
  
    
  ### DERIVATIVE ###
  
  
  ### OPTIMIZATION ###
  def optimize(self):
  
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
      pass
    elif f_format == "med": #salome
      pass
    #TODO order here the numbering ??
    return
    
  def _create_face_info_(self):
    #store the 3/4 considered face vertices and the two
    #element id which share the face
    count = 0
    temp_dict = {}
    tet_faces = [np.array([0,1,2],dtype='i2'), np.array([1,2,3],dtype='i2'),
                 np.array([0,2,3],dtype='i2'), np.array([0,1,3],dtype='i2')]
    for i,elem in enumerate(self.elements):
      elem_type = np.sum(elem != 0)
      if elem_type == 4: #Tet
        for f in tet_faces:
          key = tuple(sorted(elem[f]))
          try: temp_dict[key][1] = i
          except: 
            temp_dict[key] = [i, 0]
            count += 1
      elif elem_type == 5: #pyramid
        pass
      elif elem_type == 6: #prisms
        pass
      elif elem_type == 8: #Hex
        pass
    self.face = np.zeros((count,6),dtype='i8')
    count = 0
    for face,ids in temp_dict.items():
      if len(face) == 3:
        self.face[count,:3] = np.array(face)
        self.face[count,4] = ids[1]
        self.face[count,5] = ids[0]
      count += 1
        
    return
  
  ### INITIALIZATION FUNCTION ###
  def _initialize_all_(self):
    self._create_face_info_()
    self.cell_centers_vector = np.zeros((len(self.face),3),dtype='f8')
    self.face_normal_area = np.zeros((len(self.face),3),dtype='f8')
    self.cell_center = np.zeros((len(self.elements),3), dtype='f8')
    return
    
