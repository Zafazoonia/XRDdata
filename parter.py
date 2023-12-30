def parter(min, max):
  min_index=np.where(x==min)[0][0]
  max_index=np.where(x==max)[0][0]
  x0=x[min_index:max_index]
  y0=y[min_index:max_index]
  
  return x0,y0
