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
import matplotlib.pyplot as plt


def read_angle(f):
  res = np.genfromtxt(f)
  return res[:,2]

def plot_histograms(angle_ini, angle_final):
  #print some stats
  print("Initial stats")
  print(f"Mean non orthogonality angle: {np.mean(angle_ini)} deg")
  print(f"Max non orthogonality angle: {np.max(angle_ini)} deg")
  print("Final stats")
  print(f"Mean non orthogonality angle: {np.mean(angle_final)} deg")
  print(f"Max non orthogonality angle: {np.max(angle_final)} deg")
  #plot histogram
  fig, axarr = plt.subplots(1,2, figsize=(8,4))
  axarr[0].hist([angle_ini, angle_final], bins=20, label=("Initial", "Optimized"))
  axarr[0].set_xlabel("Non orthogonality angle")
  axarr[0].set_ylabel("Face count")
  axarr[0].legend()
  
  axarr[1].hist(angle_ini-angle_final, bins=20)
  axarr[1].set_xlabel("Delta angle")
  plt.show()
  return


if __name__ == "__main__":
  f_ini = "../src/cpp/bin/Debug/face_error_initial.txt"
  f_final = "../src/cpp/bin/Debug/face_error_final.txt"
  
  angle_ini = read_angle(f_ini)
  angle_final = read_angle(f_final)
  
  plot_histograms(angle_ini, angle_final)
  
