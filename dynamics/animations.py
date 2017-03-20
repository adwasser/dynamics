"""Swanky animations"""

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
from scipy import integrate, optimize

from .kepler import kepler_orbit

mpl.rc('animation', html='html5')

class GuidingCenterAnimation(animation.FuncAnimation):
    
    def __init__(self, a, e, nsteps=100, interval=50, **kwargs):
        """
        a : semi-major axis
        e : eccentricity
        nsteps : number of steps in theta
        kwargs : keyword arguments to pass to pyplot.subplots
        """
        self.a = a
        self.e = e
        self.nsteps = nsteps
        theta_steps = np.linspace(0, 2 * np.pi, nsteps)
        markersize = 10

        if plt.isinteractive():
            plt.ioff()
            ion = True
        # create and initialize figure and axes
        self.fig, self.ax = plt.subplots(figsize=(12, 12), **kwargs)
        self.ax.autoscale(enable=False)
        size = a * (1.3 + e)
        self.ax.set_xlim(-size, size)
        self.ax.set_ylim(-size, size)
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        
        # plot grid of guiding center and true trajectory
        theta_grid = np.linspace(0, 2 * np.pi, 100)
        r_true_grid = self.r_orbit(theta_grid, a, e)
        x_true_grid = r_true_grid * np.cos(theta_grid)
        y_true_grid = r_true_grid * np.sin(theta_grid)
        r_center_grid = self.r_orbit(theta_grid, a, 0)
        x_center_grid = r_center_grid * np.cos(theta_grid)
        y_center_grid = r_center_grid * np.sin(theta_grid)
        
        self.ax.plot(x_true_grid, y_true_grid, 'C0--', alpha=0.5)
        self.ax.plot(x_center_grid, y_center_grid, 'k-', alpha=0.5)

        plt.close()
        if ion:
            plt.ion()
            
        # save points to update later
        self.center_point = self.ax.plot(a, 0, 'ko', ms=markersize)[0]
        self.epicycle_line = self.ax.plot([], [], 'C3-.', alpha=0.5)[0]
        self.epicycle_point = self.ax.plot([], [], 'C3o', ms=markersize)[0]
        self.true_point = self.ax.plot([], [], 'C0o', ms=markersize)[0]

        super().__init__(self.fig, self.f_animate, init_func=self.f_init,
                         frames=theta_steps, interval=interval, blit=True)
    
    def rot(self, x, y, theta):
        """Rotate vector (x, y) by angle theta"""
        x_rot = x * np.cos(theta) - y * np.sin(theta)
        y_rot = x * np.sin(theta) + y * np.cos(theta)
        return x_rot, y_rot

    def r_orbit(self, theta, a, e):
        return a * (1 - e ** 2) / (1 + e * np.cos(theta))
    
    def r_epicycle(self, theta, theta_ep):
        a = self.a
        e = self.e
        a_ep = 2 * a * e
        b_ep = a * e
        x_center = a * np.cos(theta)
        y_center = a * np.sin(theta)
        x_ep = b_ep * np.cos(np.pi - theta_ep)
        y_ep = a_ep * np.sin(np.pi - theta_ep)
        x_ep, y_ep = self.rot(x_ep, y_ep, theta)
        x_ep = x_ep + x_center
        y_ep = y_ep + y_center
        return x_ep, y_ep
        
    def f_init(self):
        """Background initiation function"""
        self.center_point.set_data([], [])
        self.epicycle_line.set_data([], [])
        self.epicycle_point.set_data([], [])
        self.true_point.set_data([], [])
        return [self.center_point, self.epicycle_line, 
                self.epicycle_point, self.true_point]
    
    def f_animate(self, theta):
        """Iterative animation function"""
        a = self.a
        e = self.e
        
        # center point
        x_center = a * np.cos(theta)
        y_center = a * np.sin(theta)
        self.center_point.set_data(x_center, y_center)
        
        # epicycle line
        theta_grid = np.linspace(0, 2 * np.pi * 100)
        x_ep_grid, y_ep_grid = self.r_epicycle(theta, theta_grid)
        self.epicycle_line.set_data(x_ep_grid, y_ep_grid)
        
        # epicycle point
        x_ep, y_ep = self.r_epicycle(theta, theta)
        self.epicycle_point.set_data(x_ep, y_ep)

        # TODO true point
        x_true, y_true = kepler_orbit(theta, a, e, cart=True)
        self.true_point.set_data(x_true, y_true)
        
        return [self.center_point, self.epicycle_line, 
                self.epicycle_point, self.true_point]
