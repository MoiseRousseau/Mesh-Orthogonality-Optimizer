#!/usr/bin/env python
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


import time
from .ui import UI
import SMESH
import subprocess
import salome
import os


def export_mesh(mesh, f_out):
  out = open(f_out, 'w')
  elementCode = {4:'T', 5:'P', 6:'W', 8:'H'}
  n_nodes = mesh.NbNodes()
  n_elements = len(mesh.GetElementsByType(SMESH.VOLUME))
  out.write(f"{n_elements} {n_nodes}\n")
  for i in mesh.GetElementsByType(SMESH.VOLUME):
    nodes = mesh.GetElemNodes(i)
    out.write(elementCode[len(nodes)] + ' ')
    for x in nodes: #write
      out.write(str(x) + ' ')
    out.write('\n')
  #write node coordinates
  for i in range(1,n_nodes+1):
    X,Y,Z = mesh.GetNodeXYZ(i)
    out.write("{} {} {}\n".format(X,Y,Z))
  out.close()
  return 0 #success

def update_mesh(mesh, f_in):
  src = open(f_in, 'r')
  lines = src.readlines()
  src.close()
  for i,line in enumerate(lines):
    X, Y, Z = [float(x) for x in line.split()]
    mesh.MoveNode(i+1,X,Y,Z)
  return
  


def OrthOpt_opt(context):

  window = UI(context)
  
  result = window.show()
  
  if result:
    t = time.time()
    smesh = salome.smesh.smeshBuilder.New()
    path = '/home/%s/.config/salome/Plugins/OrthOpt/' %(os.getlogin())
  
    #get the mesh and export it
    print("\tExport mesh to OrthOpt")
    mesh = smesh.Mesh(window.mesh)
    export_mesh(mesh, path+'mesh.ugi')
    
    #get parameters
    n = window.le_penalization.text()
    ftype_code = window.menu_function.currentIndex()
    it = window.le_it.text()
    
    #optimize
    print("\tCall OrthOpt\n")
    res = subprocess.call(["./OrthOpt", "-m", "mesh.ugi", "-o", "out.xyz", "-function_type", str(ftype_code), "-penalizing_power", str(n), "-maxit", str(it)], cwd=path)
    
    #import mesh
    update_mesh(mesh, path+"out.xyz")
    if salome.sg.hasDesktop():
      salome.sg.updateObjBrowser()
      
    print ("\n\tTotal time elapsed {} s".format(time.time() - t))
      
    print ("    END \n")
    print ("####################\n\n")

  return

