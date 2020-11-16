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

import sys
from os import getcwd
sys.path.append(getcwd()+'/../src/')
from ortho_opt import Ortho_Opt
import numpy as np

import matplotlib.pyplot as plt
import time

if __name__ == "__main__":
  t = time.time()
  elements = "./10k_pts.ele"
  vertices = "./10k_pts.node"
  
  #input
  opt = Ortho_Opt()
  opt.load_elements_from_file(elements)
  opt.load_vertices_from_file(vertices)
  opt.use_numerical_derivative(False)
  opt.set_perturbation_numerical_derivative(1e-6)
  opt.set_penalizing_power(1)
  error_ini = np.copy(opt.current_face_orthogonality())
  vini = np.copy(opt.vertices)
  
  #optimize
  opt.optimize(maxiter=10)
  vfinal = np.copy(opt.vertices)
  error_final = opt.current_face_orthogonality()
  print(f"Total time: {time.time()-t}")
  
  #plot error
  plt.hist([error_ini, error_final], bins=20, label=("Initial", "Optimized"))
  plt.legend()
  plt.show()
  plt.hist(error_ini-error_final, bins=20)
  plt.show()
  
