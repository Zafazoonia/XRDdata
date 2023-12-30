#!pip install lmfit
import ast
import os
import math
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize, signal
from lmfit import models

def open_ras_file(filename):
  """Opens a .ras file and returns the data as a NumPy array."""
  with open(filename, 'r',encoding= 'unicode_escape') as f:
    data = []
    for line in f:
      if not line.startswith('*'):
        data.append([float(x) for x in line.split()])
  return np.array(data)

def open_xy_file(filename):
  """Opens a .ras file and returns the data as a NumPy array."""
  with open(filename, 'r',encoding= 'unicode_escape') as f:
    data = []
    for line in f:
      data.append([float(x) for x in line.split()])
  return np.array(data)

def parter(min, max):
  x=N[:,0]
  y=N[:,1]*N[:,2]
  min_index=np.where(x==min)[0][0]
  max_index=np.where(x==max)[0][0]
  x0=x[min_index:max_index]
  y0=y[min_index:max_index]
  return x0,y0

def update_spec_from_peaks(spec, model_indicies, peak_widths=(10, 25), **kwargs):
    x = spec['x']
    y = spec['y']
    x_range = np.max(x) - np.min(x)
    peak_indicies = signal.find_peaks_cwt(y, peak_widths)
    np.random.shuffle(peak_indicies)
    for peak_indicie, model_indicie in zip(peak_indicies.tolist(), model_indicies):
        model = spec['model'][model_indicie]
        if model['type'] in ['GaussianModel', 'LorentzianModel', 'VoigtModel']:
            params = {
                'height': y[peak_indicie],
                'sigma': x_range / len(x) * np.min(peak_widths),
                'center': x[peak_indicie]
            }
            if 'params' in model:
                model.update(params)
            else:
                model['params'] = params
        else:
            raise NotImplemented(f'model {basis_func["type"]} not implemented yet')
    print(spec)
    return peak_indicies

def generate_model(spec):
    composite_model = None
    params = None
    x = spec['x']
    y = spec['y']
    x_min = np.min(x)
    x_max = np.max(x)
    x_range = x_max - x_min
    y_max = np.max(y)
    for i, basis_func in enumerate(spec['model']):
        prefix = f'm{i}_'
        model = getattr(models, basis_func['type'])(prefix=prefix)
        if basis_func['type'] in ['GaussianModel', 'LorentzianModel', 'VoigtModel']: # for now VoigtModel has gamma constrained to sigma
            #print("ok")
            model.set_param_hint('sigma', min=1e-6, max=x_range)
            #print("sigma done")
            model.set_param_hint('center', min=x_min, max=x_max)
            #print("center done")
            model.set_param_hint('height', min=1e-6, max=1.1*y_max)
            #print("height done")
            model.set_param_hint('amplitude', min=1e-6)
            #print("amplitude done")
            # default guess is horrible!! do not use guess()
            default_params = {
                prefix+'center': x_min + x_range * random.random(),
                prefix+'height': y_max * random.random(),
                prefix+'sigma': x_range * random.random()
            }
            #print(9)
        else:
            raise NotImplemented(f'model {basis_func["type"]} not implemented yet')
            #print(8)
        if 'help' in basis_func:  # allow override of settings in parameter
            #print(7)
            for param, options in basis_func['help'].items():
                #print(6)
                model.set_param_hint(param, **options)
                #print(5)
        model_params = model.make_params(**default_params, **basis_func.get('params', {}))
        #print(4)
        if params is None:
            #print(3)
            params = model_params
            #print(2)
        else:
            params.update(model_params)
            #print(1)
        if composite_model is None:
            composite_model = model
        else:
            composite_model = composite_model + model
    return composite_model, params

def print_best_values(spec, output):
    model_params = {
        'GaussianModel':   ['amplitude', 'sigma'],
        'LorentzianModel': ['amplitude', 'sigma'],
        'VoigtModel':      ['amplitude', 'sigma', 'gamma']
    }
    best_values = output.best_values
    print('center    model   amplitude     sigma      gamma')
    for i, model in enumerate(spec['model']):
        prefix = f'm{i}_'
        values = ', '.join(f'{best_values[prefix+param]:8.3f}' for param in model_params[model["type"]])
        print(f'[{best_values[prefix+"center"]:3.3f}] {model["type"]:16}: {values}')

Filename="/home/zafazoonia/ICU/py/XRDdata/G30S10-GI-XRD-Si.ras"
N=open_ras_file(Filename)
x=N[:,0]
y=N[:,1]*N[:,2]

tablenm="/home/zafazoonia/ICU/py/XRDdata/Gete-Cubic.txt"
df=pd.read_table(tablenm, delimiter='\s+')
df.columns=['h', 'k', 'l', 'd(Å)', 'F(real)', 'F(imag)', '|F|', '2θ',
       'I', 'MID', '(λ)', 'Phase','Peak']
df['ObsInt']=df['Sigma']=df['Peak']

dfUse=df[df['2θ']<max(x)]
I=dfUse['I']
Ang=dfUse['2θ']

filename="/home/zafazoonia/ICU/py/XRDdata/Gete-Cubic.int.xy"
M=open_ras_file(filename)
Upr=np.where(M==np.ceil(max(x)))[0][0]
Lwr=np.where(M==np.ceil(min(x)))[0][0]
xx=M[Lwr:Upr,0]
yy=M[Lwr:Upr,1]*max(y)/100

fig, ax = plt.subplots()
ax.scatter(x,y,s=1)
ax.scatter(xx,yy,s=1)
pngname=os.path.splitext(Filename)[0]
plt.savefig(pngname+'.png')




for z in range(len(dfUse)):
  minx=np.ceil(Ang[z]-2)
  maxx=np.ceil(Ang[z]+2)
  x,y = parter(minx,maxx)
  fig, ax = plt.subplots()
  ax.scatter(x,y,s=1)
  #Peaks=[14.617, 24.8154, 24.4576, 28.3794, 28.9196, 29.4788, 35.8485, 38.9037,39.336, 41.3786, 44.8699, 45.9609, 47.0698]
  plt.vlines(x = Ang[z], ymin = 0, ymax = 50,
           colors = 'purple',
           label = 'vline_multiple - full height')
  spec = {
    'x': x,
    'y': y,
    'model': [
        #{'type': 'VoigtModel'},
        #{'type': 'VoigtModel'},
        #{'type': 'VoigtModel'},
        #{'type': 'VoigtModel'},
        #{'type': 'GaussianModel'},
        {'type': 'GaussianModel'},
        {'type': 'GaussianModel'},
        {'type': 'GaussianModel'},

            ]
        }

  peaks_found = update_spec_from_peaks(spec, [0, 1, 2], peak_widths=(50,))
  fig, ax = plt.subplots()
  ax.scatter(spec['x'], spec['y'], s=4)
  for i in peaks_found:
      ax.axvline(x=spec['x'][i], c='black', linestyle='dotted')
  model, params = generate_model(spec)
  output = model.fit(spec['y'], params, x=spec['x'])
  #fig, gridspec = output.plot(data_kws={'markersize':  1})
  fig, ax = plt.subplots()
  ax.scatter(spec['x'], spec['y'], s=4)
  components = output.eval_components(x=spec['x'])
  print(len(spec['model']))
  for i, model in enumerate(spec['model']):
      ax.plot(spec['x'], components[f'm{i}_'])
  model_params = {
        'GaussianModel':   ['amplitude', 'sigma'],
        'LorentzianModel': ['amplitude', 'sigma'],
        'VoigtModel':      ['amplitude', 'sigma', 'gamma']
    }
  best_values = output.best_values
  flag=0
  PkAng=x[np.where(y==max(y))[0][0]]
  for i, model in enumerate(spec['model']):
        prefix = f'm{i}_'
        values = ', '.join(f'{best_values[prefix+param]:8.3f}' for param in model_params[model["type"]])
        if flag==0:
          A=abs(PkAng-best_values[prefix+"center"])
          flag=1
        if A>abs(PkAng-best_values[prefix+"center"]):
          A=abs(PkAng-best_values[prefix+"center"])
          minI=i
  for i, model in enumerate(spec['model']):
    if i==minI:
      prefix = f'm{i}_'
      values = ', '.join(f'{best_values[prefix+param]:8.3f}' for param in model_params[model["type"]])
      dfUse['Peak'][z]=best_values[prefix+"center"]
      #print(f'[{best_values[prefix+"center"]:3.3f}]')
      numbers = [float(x) for x in values.split(',')]
      dfUse['ObsInt'][z]=numbers[0]
      dfUse['Sigma'][z]=numbers[1]

dfUse.to_csv(pngname+'.csv')      