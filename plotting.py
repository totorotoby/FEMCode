import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import sys
import time
import matplotlib.animation as animation

def plot_1d_DG(f1, f2):

    u_t = np.loadtxt(f1)
    x = np.loadtxt(f2)

    to_plot = []
    for i in range(len(u_t)):
        u = u_t[i]
        line = []
        for j in range(0, len(u)-1, 2):
            line.append([(x[j], u[j]), (x[j+1], u[j+1])])        
        to_plot.append(line)
            
    return to_plot, u_t


def animate(i):


    min_u = np.amin(u_t)
    max_u = np.amax(u_t)
    diff = max_u - min_u
    
    segs = LineCollection(to_plot[i])
    ax.clear()
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(min_u, max_u + (.1*diff))
    ax.add_collection(segs)
    plt.show()
    

to_plot, u_t = plot_1d_DG(sys.argv[1], sys.argv[2]) 

fig, ax = plt.subplots()

animation.FuncAnimation(fig, animate, frames=range(1, len(u_t)), interval=1000, repeat=False)

plt.show()


