#task 6

import numpy as np

""""" 
    Calculating the refractive index using snells law. 
    Args - input parameters of function. 
    direc - direction vector, normal - of surface, n_1,n_2 - refecative indexes.
    
    Returns the refracted direction vecotor or none if total internal reflection.
    """

def refract(direc, normal, n_1, n_2):
    
    # Normalize the direction and normal vectors
    direc = np.array(direc) / np.linalg.norm(direc)
    normal = np.array(normal) / np.linalg.norm(normal)

    # Calculate the cosine of the angle of incidence using dot product
    cos_theta_i = -np.dot(direc, normal)

    # Calculate the sine of the angle of incidence using Pythagorean identity
    sin_theta_i_squared = 1 - cos_theta_i**2
    sin_theta_t_squared = (n_1 / n_2)**2 * sin_theta_i_squared

    # Check for total internal reflection
    if sin_theta_t_squared > 1:
        return None

    # Calculate the cosine of the angle of refraction using Pythagorean identity
    cos_theta_t = np.sqrt(1 - sin_theta_t_squared)

    # Calculate the refracted direction
    refracted_direc_parallel = (n_1 / n_2) * (direc + cos_theta_i * normal)
    refracted_direc_perpendicular = -cos_theta_t * normal
    refracted_direc = refracted_direc_parallel + refracted_direc_perpendicular

    # Normalize the refracted direction to ensure it is a unit vector
    refracted_direc = refracted_direc / np.linalg.norm(refracted_direc)

    return refracted_direc



