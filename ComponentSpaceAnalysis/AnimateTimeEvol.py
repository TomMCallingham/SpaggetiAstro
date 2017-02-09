import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from ComponentSpaceAnalysis.TimeEvolution import Evolution
from numba import jit



fig=plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(0, 2*np.pi), ylim=(-2, 2))
line, = ax.plot([],[])
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def animate(t):
    #currentx=X
    #currenty=Y[t,:]


    line.set_data(X,Y[t,:])
    time_text.set_text(time_template % t)
    return line,time_text



ani=animation.FuncAnimation(fig,animate,frames=100,interval=20,init_func=init,blit=True)
plt.show()
