import numpy as np
import os
import yaml
import torch
np.arange(10)

os.getcwd()

#!pip install matplotlib
#!pip install torch
#pip install neuralhydrology
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import torch

import netCDF4

from neuralhydrology.evaluation import metrics
from neuralhydrology.nh_run import start_run, eval_run

# by default we assume that you have at least one CUDA-capable NVIDIA GPU or MacOS with Metal support
if torch.cuda.is_available() or torch.backends.mps.is_available():
    start_run(config_file=Path("allbasins-ens3.yml"))

# fall back to CPU-only mode
else:
    start_run(config_file=Path("allbasins-ens3.yml"), gpu=-1)