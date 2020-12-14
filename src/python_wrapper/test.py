
from OrthOpt import OrthOpt
import numpy as np

if __name__ == "__main__":
  f_v = "../../examples/10k_pts.node"
  f_e = "../../examples/10k_pts.ele"
  vertices = np.genfromtxt(f_v)
  elements = np.genfromtxt(f_e)
  
  opt = OrthOpt(vertices, elements)
