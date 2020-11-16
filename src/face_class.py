
import numpy as np

class Face:
  def __init__(self, unique_id=-999):
    self.unique_id = unique_id
    self.element_id_up = -999 #element id in the direction of face normal
    self.element_id_dn = -999
    self.vertices = -999 #the vertices of the face
    self.vertice_up = -999 #the F position vertice
    self.vertice_dn = -999 #the E position vertice
    self.type = -999
    self.area = -999.
    self.normal = np.zeros(3,dtype='f8')
    self.cell_centers_vector = np.zeros(3,dtype='f8')
    self.cell_centers_vector_norm = -999.
    self.error = -999.
    self.weight = 1.
    return
    
  def populate_face(self, id_up, id_dn, vertices, vertice_up, vertice_dn):
    self.element_id_up = id_up
    self.element_id_dn = id_dn
    self.vertices = vertices
    self.vertice_up = vertice_up
    self.vertice_dn = vertice_dn
    self.type = len(vertices)
    return
    
  def print_informations(self):
    print(f"Element id up: {self.element_id_up}")
    print(f"Element id dn: {self.element_id_dn}")
    print("Vertices:", end='')
    for x in self.vertices: print(f" {x},", end='')
    print(f"\nVertice up: {self.vertice_up}")
    print(f"Vertice dn: {self.vertice_dn}")
    if self.type == 3: print("Face type: triangular")
    else: print("Face type: quad")
    print(f"Area: {self.area}")
    print(f"Cell center vector across face: {self.cell_centers_vector}")
    print(f"Cell center vector across face norm: {self.cell_centers_vector_norm}")
    print(f"Normal: {self.normal}\n")
    return
    
  def compute_error(self):
    self.error = 1 - np.dot(self.cell_centers_vector, self.normal)
    return self.error


class Vertice:
  def __init__(self, coor, index=-999):
    self.coor = coor
    self.index = index
    self.implied_in_faces = {}
    return
    
  def add_face_and_position(self, face_id, position):
    #position could be A, E or F
    self.implied_in_faces[face_id] = position
    return
    
