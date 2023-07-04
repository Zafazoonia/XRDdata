import numpy as np
import matplotlib.pyplot as plt

def open_ras_file(filename):
  """Opens a .ras file and returns the data as a NumPy array."""
  with open(filename, 'r',encoding= 'unicode_escape') as f:
    data = []
    for line in f:
      if not line.startswith('*'):
        data.append([float(x) for x in line.split()])
  return np.array(data)

def RasRead(filename):
  N=open_ras_file(filename)
  x=N[:,0]
  y=N[:,1]*N[:,2]
  fig, ax = plt.subplots()
  ax.scatter(x,y,s=1)
  return x,y
