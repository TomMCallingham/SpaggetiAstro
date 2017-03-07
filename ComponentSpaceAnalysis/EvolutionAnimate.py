import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from ComponentSpaceAnalysis.Evolution import *
#savename='DistSaveM7T6t4J4N200.npy'
savename='DistTempSave.npy'
Dist = np.load(savename)
N = np.size(Dist[0, :, 0])
ShowT = np.size(Dist[0, 0, :])
H = np.linspace(0, 1, N)

#TotalVol=1e10


fig=plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-0.1, 1.1), ylim=(0, 2 * TotalVol))
plt.ylabel('Total Vol')
plt.xlabel('H, Line position')
line, = ax.plot([],[])
time_template = 'time = %s e4years'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def animate(t):

    line.set_data(H,Dist[0,:,t])
    time_text.set_text(time_template % t)
    return line,time_text



anim=animation.FuncAnimation(fig,animate,frames=ShowT,interval=100,init_func=init,blit=True)
anim.save('Anim%s.mp4'%savename)#fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()
