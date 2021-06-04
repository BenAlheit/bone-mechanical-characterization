import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib

SMALL_SIZE = 8
MEDIUM_SIZE = 11
BIGGER_SIZE = 13
FONT = 16

plt.rc('font', size=FONT)          # controls default text sizes
plt.rc('axes', titlesize=FONT)     # fontsize of the axes title
plt.rc('axes', labelsize=FONT)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}',
    r'\usepackage{amssymb}']

