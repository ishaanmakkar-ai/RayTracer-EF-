"""
This module contains the Ray class for simulating optical rays.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Ray:
    """
    Class to represent an optical ray with multiple points and directions in 3D space.
    
    The Ray class keeps track of the path of a ray through space by storing its 
    positions and directions as it propagates.
    """

    def __init__(self, pos = None, direc = None):
        """
        Initialize the ray with a starting position and direction.
        
        Args:
            pos (list or np.ndarray, optional): Starting position of the ray. 
                                                Defaults to [0, 0, 0].
            direc (list or np.ndarray, optional): Starting direction of the ray. 
                                                  Defaults to [0, 0, 1].
        """
        if pos is None:
            pos =[0., 0., 0.]
        if direc is None:
            direc =[0., 0., 1.]

        self._pos = np.array(pos)
        self._direc = np.array(direc)
        self._points = [self._pos]

 
    def pos(self):
        """
        Return the current position of the ray.
        
        Returns:
            np.ndarray: The current position of the ray.
        """
        return self._pos

  
    def direc(self):
        """
        Return the current direction of the ray.
        
        Returns:
            np.ndarray: The current direction of the ray.
        """
        return self._direc

    def append(self, pos, direc):
        """
        Append a new position and direction to the ray.
        
        Args:
            pos (list or np.ndarray): New position to add to the ray.
            direc (list or np.ndarray): New direction associated with the new position.
        """
        self._pos = np.array(pos)
        self._direc = np.array(direc)
        self._points.append(self._pos)

  
    def vertices(self):
        """
        Return all positions along the ray.
        
        Returns:
            list of np.ndarray: All positions the ray has passed through.
        """
        return self._points

class RayBundle:
    def __init__(self, rmax=5., nrings=5, multi=6):
        self.rays = []
        for r in np.linspace(0, rmax, nrings):
            for theta in np.linspace(0, 2 * np.pi, multi, endpoint=False):
                x, y = r * np.cos(theta), r * np.sin(theta)
                self.rays.append(Ray(pos=[x, y, 0], direc=[0, 0, 1]))

    def propagate_bundle(self, elements):
        for ray in self.rays:
            for elem in elements:
                elem.propagate_ray(ray)

    def track_plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')#check this not sure what it means
        for ray in self.rays:
            vertices = np.array(ray.vertices())
            ax.plot(vertices[:, 0], vertices[:, 1], vertices[:, 2])
        return fig
    
    def rms(self):
        positions = np.array([ray.pos() for ray in self.rays])
        rms_value = np.sqrt(np.mean(np.sum(positions[:, :2]**2, axis=1)))
        return rms_value


    def spot_plot(self):
        fig, ax = plt.subplots()
        positions = np.array([ray.pos() for ray in self.rays])
        ax.plot(positions[:, 0], positions[:, 1], 'o')
        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')
        ax.set_title('Spot Diagram')
        ax.grid(True)
        return fig


