import os
import matplotlib.pyplot as plt
import numpy as np
import time

# Global figure variable
# This is to make sure each plot is drawn in a new window, no matter which plotting methods are used
n_fig = 1

def basic_xy(x,y,color='b'):
    
    global n_fig
    figure = plt.figure(n_fig)
    figure.add_subplot(1, 1, 1, axisbg='1') # Change background color here
    plt.gca().set_aspect('equal')  
    plt.plot(x,y,color)
    plt.show()
    
#    time.sleep(5)
    n_fig += 1

def body_wake_plot(Swimmers):
    
    global n_fig
    figure = plt.figure(n_fig)
    plt.clf()
    figure.add_subplot(1, 1, 1, axisbg='1') # Change background color here
    plt.gca().set_aspect('equal')
    maxpercentile = 95 # For truncating outliers
    
    # Gather circulations of all swimmers into a color array
    color = []
    for Swim in Swimmers:
        Swim.n_color = len(Swim.Wake.gamma[1:-1])
        color = np.append(color, Swim.Wake.gamma[1:-1])
        Swim.i_color = len(color)-Swim.n_color
    
    # Make color map based on vorticity
    # Take a look at positive and negative circulations separately
    if np.min(color) < 0: # Check if negative circulations exist (in case of short simulations)
        # Truncate any negative outliers
        color[color < np.percentile(color[color < 0], 100-maxpercentile)] = np.percentile(color[color < 0], 100-maxpercentile)
        # Normalize negative circulations to [-1,0)
        color[color < 0] = -color[color < 0]/np.min(color)
    if np.max(color) > 0: # Check if positive circulations exist (in case of short simulations)
        # Truncate any positive outliers
        color[color > np.percentile(color[color > 0], maxpercentile)] = np.percentile(color[color > 0], maxpercentile)
        # Normalize positive circulations to (0,1]
        color[color > 0] = color[color > 0]/np.max(color)
        
    for Swim in Swimmers:
        # Extract color map for the individual Swim
        c = color[Swim.i_color:Swim.i_color+Swim.n_color]
        # Scatter plot of wake points with red-white-blue colormap, as well as body outline and edge panel segment
        plt.scatter(Swim.Wake.x[1:-1], Swim.Wake.z[1:-1], s=30, c=c, edgecolors='none', cmap=plt.get_cmap('bwr_r'))
        plt.plot(Swim.Body.AF.x, Swim.Body.AF.z, 'k')
        plt.plot(Swim.Edge.x, Swim.Edge.z, 'g')
        plt.show()
    
    n_fig += 1
    
def cp_plot(Body):
    
    global n_fig
    figure = plt.figure(n_fig)
    figure.add_subplot(1, 1, 1, axisbg='1') # Change background color here
    plt.gca().set_aspect('equal')
    plt.gca().invert_yaxis()
    
    plt.plot(Body.AF.x_col[:Body.N/2], Body.cp[:Body.N/2], 'g')
    plt.plot(Body.AF.x_col[Body.N/2:], Body.cp[Body.N/2:], 'b')
    plt.plot(Body.AF.x, -Body.AF.z, 'k')
    plt.plot(Body.AF.x_col, -Body.AF.z_col, 'r')
    
    n_fig += 1
    
def drag_vs_period(Body,RHO,t):
    
    global n_fig
    figure = plt.figure(n_fig)
    figure.add_subplot(1, 1, 1, axisbg='1') # Change background color here
    plt.xlabel('tau')
    plt.ylabel('Coefficent of drag')
    
    plt.plot(t[4:]*Body.F, -Body.drag[3:]/(0.5*RHO*Body.V0**2), 'b')
    
    n_fig += 1
    
def lift_vs_period(Body,RHO,t):
    
    global n_fig
    figure = plt.figure(n_fig)
    figure.add_subplot(1, 1, 1, axisbg='1') # Change background color here
    plt.xlabel('tau')
    plt.ylabel('Coefficent of lift')
    
    plt.plot(t[4:]*Body.F, -Body.lift[3:]/(0.5*RHO*Body.V0**2), 'g')
    
    n_fig += 1
    
def plot_n_go(Edge, Body, Solid, V0, T, HEAVE):
    global n_fig
    # Determine if the output directory exists. If not, create the directory.
    if not os.path.exists('./movies'):
        os.makedirs('./movies')
        
    figure = plt.figure(1)
    figure.add_subplot(1, 1, 1, axisbg='1') # Change background color here
#    plt.gca().set_aspect('equal')
    plt.gca().invert_yaxis()
    
#    plt.plot(Body.AF.x_col[:Body.N/2], Body.cp[:Body.N/2]/100, 'g')
#    plt.plot(Body.AF.x_col[Body.N/2:], Body.cp[Body.N/2:]/100, 'b')
    plt.plot(Body.AF.x, Body.AF.z, 'k')
#    plt.plot(Body.AF.x_col, Body.AF.z_col, 'r')
    
#    plt.plot(Solid.tempNodes[:,0], Solid.tempNodes[:,1], 'k')
    
    plt.xlim((np.min(Body.AF.x)-0.02, np.min(Body.AF.x)+0.22))
#    plt.plot(Edge.x, Edge.z, 'g')
    plt.ylim((-0.05, 0.05))
    
    rigidX = np.zeros(2)
    rigidZ = np.zeros(2)
    rigidX[0] = 0.
    rigidX[1] = 0.2
    rigidZ[0] = 0.
    rigidZ[1] = 0.
    rigidX += V0*T
    rigidZ += HEAVE
    plt.plot(rigidX, rigidZ, 'm')
    
    
    figure.savefig('./movies/%05i.png' % (n_fig), format='png')
    plt.clf()
    
    n_fig += 1  

def body_plot(Edge, Body):
    global n_fig
    
    # Determine if the output directory exists. If not, create the directory.
    if not os.path.exists('./movies'):
        os.makedirs('./movies')
        
    figure = plt.figure(1)
    figure.add_subplot(1, 1, 1, axisbg='1') # Change background color here
#    plt.gca().set_aspect('equal')
    plt.gca().invert_yaxis()
    
#    plt.plot(Body.AF.x_col[:Body.N/2], Body.cp[:Body.N/2]/100, 'g')
#    plt.plot(Body.AF.x_col[Body.N/2:], Body.cp[Body.N/2:]/100, 'b')
    plt.plot(Body.AF.x, Body.AF.z, 'k')
    plt.plot(Body.AF.x_col, Body.AF.z_col, 'r')

    plt.xlim((np.min(Body.AF.x)-0.02, np.min(Body.AF.x)+0.22))
    plt.plot(Edge.x, Edge.z, 'g')
    plt.ylim((-0.05, 0.05))
    
    figure.savefig('./movies/%05i.png' % (n_fig), format='png')
    plt.clf()
    
    n_fig += 1

def body(x,y,color='b'):
    figure = plt.figure(2)
    figure.add_subplot(1, 1, 1, axisbg='1') # Change background color here
    plt.gca().set_aspect('equal')  
    plt.plot(x,y,color)
    plt.xlim((np.min(x)-0.02, np.min(x)+0.22))
    plt.ylim((-0.05, 0.05))
    plt.show()