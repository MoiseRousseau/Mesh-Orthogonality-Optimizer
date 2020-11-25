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
from OrthOpt import Ortho_Opt
import numpy as np

import matplotlib.pyplot as plt
import time


def optimize(f_vertices, f_elements):
  t = time.time()
  
  #input
  opt = Ortho_Opt()
  opt.load_elements_from_file(f_elements)
  opt.load_vertices_from_file(f_vertices)
  opt.use_numerical_derivative(False)
  opt.set_perturbation_numerical_derivative(1e-6)
  opt.set_penalizing_power(1)
  error_ini = np.arccos(np.abs(1-opt.current_face_orthogonality())) * 180/np.pi
  
  #optimize
  opt.optimize(maxiter=10)
  vfinal = np.copy(opt.vertices)
  error_final = np.arccos(np.abs(1-opt.current_face_orthogonality())) * 180/np.pi
  print(f"Total time: {time.time()-t}")
  
  #plot error
  fig, axarr = plt.subplots(1,2, figsize=(8,4))
  axarr[0].hist([error_ini, error_final], bins=20, label=("Initial", "Optimized"))
  axarr[0].set_xlabel("Non orthogonality angle")
  axarr[0].set_ylabel("Face count")
  axarr[0].legend()
  
  axarr[1].hist(error_ini-error_final, bins=20)
  axarr[1].set_xlabel("Delta angle")
  plt.show()



if __name__ == "__main__":
  #TODO command line argument
  exit(0)
  
