import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from functions_general import point_vectors, panel_vectors
import numpy as np

# Global figure variable
# This is to make sure each plot is drawn in a new window, no matter which plotting methods are used
n_fig = 1

def basic_xyz(x,y,z,color='b'):
    
    global n_fig
    fig = plt.figure(n_fig)
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')
    ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-0.2, 2.2)
    ax.set_zlim(-1.2, 1.2)
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
    
def plot_n_go(Swimmers, V0, T, HEAVE, i, SW_PLOT_FIG):
    global n_fig
    
    if SW_PLOT_FIG:
        figure = plt.figure(1)
        figure.add_subplot(1, 1, 1, axisbg='1') # Change background color here
        figure.set_size_inches(16, 9)
        plt.gca().set_aspect('equal')
        plt.tick_params(labelsize=28)
        plt.xticks(np.arange(-15.0, 15.0, 0.2))
        maxpercentile = 95 # For truncating outliers
        
        if (i > 1):
            # Gather circulations of all swimmers into a color array
            color = []
            for Swim in Swimmers:
                Swim.n_color = len(Swim.Wake.gamma[1:i])
                color = np.append(color, Swim.Wake.gamma[1:i])
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
            if (i > 1):
                # Extract color map for the individual Swim
                c = color[Swim.i_color:Swim.i_color+Swim.n_color]
                # Scatter plot of wake points with red-white-blue colormap, as well as body outline and edge panel segment
    #            for idx in xrange(i):
                plt.scatter(Swim.Wake.x[1:i], Swim.Wake.z[1:i], s=30, c=c, edgecolors='none', cmap=plt.get_cmap('bwr_r'))
    #            plt.scatter(Swim.Wake.x[1:-1], Swim.Wake.z[1:-1], s=30, c=c, edgecolors='none', cmap=plt.get_cmap('bwr_r'))
            plt.plot(Swim.Body.AF.x, Swim.Body.AF.z, 'k')
            plt.plot(Swim.Edge.x, Swim.Edge.z, 'g')
    
        # Determine if the output directory exists. If not, create the directory.
        if not os.path.exists('./movies'):
            os.makedirs('./movies')
        
        plt.axis([np.min(Swim.Body.AF.x)-0.05, np.min(Swim.Body.AF.x)+1.75, -0.5, 0.5])
#        plt.axis([np.min(Swim.Body.AF.x)-0.75, np.min(Swim.Body.AF.x)+25.5, -7.5, 7.5])
        plt.xlabel('$X$ $[m]$', fontsize=28)
        plt.ylabel('$Z$ $[m]$', fontsize=28)
        
        plt.axes([0.125, 0.677, 0.2, 0.2])
        plt.gca().set_aspect('equal')
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
        plt.axis([np.min(Swim.Body.AF.x)-0.05, np.min(Swim.Body.AF.x)+0.15, -0.06, 0.06])
#        plt.axis([np.min(Swim.Body.AF.x)-0.75, np.min(Swim.Body.AF.x)+2.25, -0.9, 0.9])
        for Swim in Swimmers:
            plt.plot(Swim.Body.AF.x, Swim.Body.AF.z, 'k')
            plt.plot(Swim.Edge.x, Swim.Edge.z, 'g')
            
        figure.savefig('./movies/%05i.png' % (n_fig), format='png')
        
        plt.clf()
    
    n_fig += 1  

def body_plot(Swimmers, SW_PLOT_FIG):
    global n_fig
        
    if SW_PLOT_FIG:
        figure = plt.figure(1)
        ax = figure.gca(projection='3d')
#        ax.set_aspect('equal') This feature has not been implemented in 3D plotting yet
        figure.set_size_inches(16, 9)
        plt.tick_params(labelsize=28)
            
        for Swim in Swimmers:
            Nc = Swim.Body.BF.x.shape[0]
            Ns = Swim.Body.BF.x.shape[1]
            (nx, ny, nz, txs, tys, tzs, txc, tyc, tzc, lps, lpc) = panel_vectors(Swim.Body.AF.x, Swim.Body.AF.y, Swim.Body.AF.z, Ns, Nc)
            x = Swim.Body.AF.x_mid[::5,::5] + 0.1 * nx[::5,::5]
            y = Swim.Body.AF.y_mid[::5,::5] + 0.1 * ny[::5,::5]
            z = Swim.Body.AF.z_mid[::5,::5] + 0.1 * nz[::5,::5]
            nx = nx[::5,::5]
            ny = ny[::5,::5]
            nz = nz[::5,::5]
            ax.quiver(x, y, z, nx, ny, nz, length=0.1, color='r')
            ax.quiver(Swim.Body.AF.x_mid[::5,::5], Swim.Body.AF.y_mid[::5,::5], Swim.Body.AF.z_mid[::5,::5], txs[::5,::5], tys[::5,::5], tzs[::5,::5], length=0.1, color='g')
            ax.quiver(Swim.Body.AF.x_mid[::5,::5], Swim.Body.AF.y_mid[::5,::5], Swim.Body.AF.z_mid[::5,::5], txc[::5,::5], tyc[::5,::5], tzc[::5,::5], length=0.1, color='b')
            ax.plot_surface(Swim.Body.AF.x, Swim.Body.AF.y, Swim.Body.AF.z, rstride=1, cstride=1, linewidth=0, color='k', antialiased=True)
    
        # Determine if the output directory exists. If not, create the directory.
        if not os.path.exists('./movies'):
            os.makedirs('./movies')
        
        ax.set_xlim(np.min(Swim.Body.AF.x)-1.25, np.min(Swim.Body.AF.x)+1.25)
        ax.set_ylim(-1.25, 1.25)
        ax.set_zlim(-1.25, 1.25)
        ax.view_init(elev=18, azim=-124)

        ax.set_xlabel('$X$ $[m]$', fontsize=28)
        ax.set_ylabel('$Y$ $[m]$', fontsize=28)
        ax.set_zlabel('$Z$ $[m]$', fontsize=28)
        ax.xaxis._axinfo['label']['space_factor'] = 2.0
        ax.yaxis._axinfo['label']['space_factor'] = 2.0
        ax.zaxis._axinfo['label']['space_factor'] = 2.0
            
        figure.savefig('./movies/%05i.png' % (n_fig), format='png')
        
        plt.clf()
    
    n_fig += 1
    
def body_wake_plot(Swimmers, SW_PLOT_FIG, i):
    global n_fig

    if SW_PLOT_FIG:
        figure = plt.figure(1)
        ax = figure.gca(projection='3d')
#        ax.set_aspect('equal') This feature has not been implemented in 3D plotting yet
        figure.set_size_inches(16, 9)
        plt.tick_params(labelsize=28)
            
        for Swim in Swimmers:
            ax.plot_surface(Swim.Body.AF.x, Swim.Body.AF.y, Swim.Body.AF.z, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False)
            if (i > 1):
                ax.plot_surface(Swim.Wake.x[1:i], Swim.Wake.y[1:i], Swim.Wake.z[1:i], rstride=1, cstride=1, alpha=0.3)
    
        # Determine if the output directory exists. If not, create the directory.
        if not os.path.exists('./movies'):
            os.makedirs('./movies')
        
        ax.set_xlim(-4.25, 1.25)
        ax.set_ylim(-2.75, 2.75)
        ax.set_zlim(-2.75, 2.75)
        ax.view_init(elev=18, azim=-124)

        ax.set_xlabel('$X$ $[m]$', fontsize=28)
        ax.set_ylabel('$Y$ $[m]$', fontsize=28)
        ax.set_zlabel('$Z$ $[m]$', fontsize=28)
        ax.xaxis._axinfo['label']['space_factor'] = 2.0
        ax.yaxis._axinfo['label']['space_factor'] = 2.0
        ax.zaxis._axinfo['label']['space_factor'] = 2.0
            
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