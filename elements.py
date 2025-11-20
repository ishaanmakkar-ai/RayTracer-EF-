
#Task 4
import numpy as np
from .physics import refract #task 7


class OpticalElement:
    def intercept(self, ray):
        raise NotImplementedError('intercept() needs to be implemented in derived classes')

    def propagate_ray(self, ray):
        raise NotImplementedError('propagate_ray() needs to be implemented in derived classes')

class SphericalRefraction(OpticalElement):
    def __init__(self, z_0=10, aperture=5., curvature=0.02, n_1=1., n_2=1.5):
        self._z_0 = z_0
        self._aperture = aperture
        self._curvature = curvature
        self._n_1 = n_1
        self._n_2 = n_2
        self._radius = 1 / curvature if curvature != 0 else np.inf


    def z_0(self):
        return self._z_0

    def aperture(self):
        return self._aperture

    
    def curvature(self):
        return self._curvature

    
    def n_1(self):
        return self._n_1

    
    def n_2(self):
        return self._n_2
    
    #task 11
    def focal_point(self):
        f = self._n_2 / ((self._n_2 - self._n_1) * self._curvature)
        return self._z_0 + f
    

#Task 5

    def intercept(self, ray):
        ray_origin = ray.pos()
        ray_direction = ray.direc()

        if self._curvature == 0:
            # Planar surface
            t = (self._z_0 - ray_origin[2]) / ray_direction[2]
            if t < 0:
                return None
            intersection = ray_origin + t * ray_direction
            if np.linalg.norm(intersection[:2]) <= self._aperture:
                return intersection
            else:
                return None

        # Center of the spherical surface
        sphere_center = np.array([0, 0, self._z_0 + self._radius])

        # Vector from sphere center to ray origin
        oc = ray_origin - sphere_center

        # Coefficients for the quadratic equation
        b = 2 * np.dot(oc, ray_direction)
        c = np.dot(oc, oc) - self._radius**2

        # Discriminant of the quadratic equation
        discriminant = b**2 - 4 * c

        if discriminant < 0:
            return None  # No intersection

        # Calculate the two possible intersections (t1 and t2)
        sqrt_discriminant = np.sqrt(discriminant)
        t1 = (-b - sqrt_discriminant) / 2
        t2 = (-b + sqrt_discriminant) / 2

        # Choose the smallest positive intersection
        if t1 > 0 and (t1 < t2 or t2 <= 0):
            t = t1
        elif t2 > 0:
            t = t2
        else:
            return None  # Both intersections are negative

        # Calculate the intersection point
        intersection = ray_origin + t * ray_direction

        # Check if the intersection point is within the aperture
        if np.linalg.norm(intersection[:2]) <= self._aperture:
            return intersection
        else:
            return None

# Testing the SphericalRefraction class - may not be needed. check when redoing
# moved out of class

#Task 7 


    def propagate_ray(self, ray):
        """
        Ray propagating through new surface, appending new position and direction to ray
        """
        intersection = self.intercept(ray)
        if intersection is None:
            return #None #no valid intersection
        normal = intersection - np.array([0,0, self._z_0 + self._radius])
        normal = normal / np.linalg.norm(normal)#normalise

        new_direction = refract(ray.direc(), normal, self._n_1, self._n_2)
        if new_direction is None:
            return #None #total internal reflection. don't update the ray.
    
        ray.append(intersection, new_direction)
        
        
#task 9        
     
class OutputPlane(OpticalElement):
    def __init__(self, z_0):
        self._z_0 = z_0

    def z_0(self):
        return self._z_0
    
    def intercept(self, ray):
        pos = ray.pos()
        direc = ray.direc()

        if direc[2] == 0:
            return None
        
        t = (self._z_0 - pos[2])/direc[2]
        if t < 0:
            return None
        intercept_point = pos + t * direc
        return intercept_point
    
    def propagate_ray(self, ray):
        intercept_point = self.intercept(ray)
        if intercept_point is not None:
            ray.append(intercept_point, ray.direc())









