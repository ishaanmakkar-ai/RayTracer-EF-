"""Analysis module."""
import matplotlib.pyplot as plt
import numpy as np
from raytracer.elements import SphericalRefraction
from raytracer.elements import OutputPlane
from raytracer.rays import Ray, RayBundle



def task8():
    """
    Task 8.

    In this function you should check your propagate_ray function properly
    finds the correct intercept and correctly refracts a ray. Don't forget
    to check that the correct values are appended to your Ray object.
    """
    #range of initial rays
    
    ray1 = Ray()
    ray2 = Ray(pos = [1,0,-10])
    #ray3 = Ray(pos = [-1,0,-10], direc = [0,0,1]),
    #ray4 = Ray(pos = [0,1,-10], direc = [0,0,1]),
    #ay5 = Ray(pos = [0,-1,-10], direc = [0,0,1])
    

    #spherical refraction element
    sr = SphericalRefraction( z_0=10, aperture=5, curvature= 0.02, n_1=1.0, n_2= 1.5)
    sr.propagate_ray(ray1)
    sr.propagate_ray(ray2)
    #print results if needed?

def task10():#check task 10 no plot made
    """
    Task 10.

    In this function you should create Ray objects with the given initial positions.
    These rays should be propagated through the surface, up to the output plane.
    You should then plot the tracks of these rays.
    This function should return the matplotlib figure of the ray paths.

    Returns:
        Figure: the ray path plot.
    """
    #defining optical elememts - given in the task
    spherical_surface = SphericalRefraction(z_0=100, curvature = 0.03, n_1= 1.0, n_2= 1.5, aperture=10)
    output_plane = OutputPlane(z_0 = 250)
#mind the commas
    initial_positions = [
        [0, 4, 0],
        [0, 1, 0],
        [0, 0.2, 0],
        [0, 0, 0],
        [0, -0.2, 0],
        [0, -1, 0],
        [0, -4, 0],
    ]

    #check
    fig, ax = plt.subplots()
    for pos in initial_positions:
        ray = Ray(pos = pos, direc = [0, 0, 1])
        elemenents = [spherical_surface, output_plane]
        for elem in elemenents:
            elem.propagate_ray(ray)


        vertices = ray.vertices()
        vertices = np.array(vertices)
        ax.plot(vertices[:, 2], vertices[:, 1]) #, label = f"Initial y={pos[1]}")

    ax.set_xlabel('z (mm)')
    ax.set_ylabel('y (mm)')
    #ax.title
    ax.legend()
    ax.grid(True)

    return fig


def task11():
    """
    Task 11.

    In this function you should propagate the three given paraxial rays through the system
    to the output plane and the tracks of these rays should then be plotted.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for ray paths
    2. the calculated focal point.

    Returns:
        tuple[Figure, float]: the ray path plot and the focal point
    """
    #spherical surface with parameters
    spherical_surface = SphericalRefraction(z_0=100, curvature = 0.03, n_1= 1.0, n_2= 1.5, aperture=10)
    #focal point
    focal_point = spherical_surface.focal_point()

    output_plane = OutputPlane(z_0=focal_point)

    #initial positons
    initial_positions = [
        [0.1,0.1,0.],
        [0.,0.,0.],
        [-0.1,-0.1,0.],
    ]
    fig, ax = plt.subplots()
    for pos in initial_positions:
        ray = Ray(pos = pos, direc = [0, 0, 1])
        elemenents = [spherical_surface, output_plane]
        for elem in elemenents:
            elem.propagate_ray(ray)


        vertices = ray.vertices()
        vertices = np.array(vertices)
        ax.plot(vertices[:, 2], vertices[:, 1]) #, label = f"Initial y={pos[1]}")

    ax.set_xlabel('z (mm)')
    ax.set_ylabel('y (mm)')
    #ax.title
    ax.legend()
    ax.grid(True)

    return fig, focal_point

    


def task12():
    """
    Task 12.

    In this function you should create a RayBunble and propagate it to the output plane
    before plotting the tracks of the rays.
    This function should return the matplotlib figure of the track plot.

    Returns:
        Figure: the track plot.

    """
    spherical_surface = SphericalRefraction(z_0=100, curvature=0.03, n_1=1.0, n_2=1.5, aperture=10)
    focal_point = spherical_surface.focal_point()
    output_plane = OutputPlane(z_0=focal_point)
    
    ray_bundle = RayBundle()
    
    elements = [spherical_surface, output_plane]
    ray_bundle.propagate_bundle(elements)
    
    fig = ray_bundle.track_plot()
    return fig



def task13():
    """
    Task 13.

    In this function you should again create and propagate a RayBundle to the output plane
    before plotting the spot plot.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the spot plot
    2. the simulation RMS

    Returns:
        tuple[Figure, float]: the spot plot and rms
    """

    spherical_surface = SphericalRefraction(z_0=100, curvature=0.03, n_1=1.0, n_2=1.5, aperture=10)
    focal_point = spherical_surface.focal_point()
    output_plane = OutputPlane(z_0=focal_point)

    ray_bundle = RayBundle()
    elements = [spherical_surface, output_plane]
    ray_bundle.propagate_bundle(elements)

    fig = ray_bundle.spot_plot()
    rms_value = ray_bundle.rms()
    
    return fig, rms_value


def task14():
    """
    Task 14.

    In this function you will trace a number of RayBundles through the optical system and
    plot the RMS and diffraction scale dependence on input beam radii.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the diffraction scale plot
    2. the simulation RMS for input beam radius 2.5
    3. the diffraction scale for input beam radius 2.5

    Returns:
        tuple[Figure, float, float]: the plot, the simulation RMS value, the diffraction scale.
    """
    return


def task15():
    """
    Task 15.

    In this function you will create plano-convex lenses in each orientation and propagate a RayBundle
    through each to their respective focal point. You should then plot the spot plot for each orientation.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the spot plot for the plano-convex system
    2. the focal point for the plano-convex lens
    3. the matplotlib figure object for the spot plot for the convex-plano system
    4  the focal point for the convex-plano lens


    Returns:
        tuple[Figure, float, Figure, float]: the spot plots and rms for plano-convex and convex-plano.
    """
    return


def task16():
    """
    Task 16.

    In this function you will be again plotting the radial dependence of the RMS and diffraction values
    for each orientation of your lens.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the diffraction scale plot
    2. the RMS for input beam radius 3.5 for the plano-convex system
    3. the RMS for input beam radius 3.5 for the convex-plano system
    4  the diffraction scale for input beam radius 3.5

    Returns:
        tuple[Figure, float, float, float]: the plot, RMS for plano-convex, RMS for convex-plano, diffraction scale.
    """
    return


if __name__ == "__main__":

    # Run task 8 function
    task8()

    # Run task 10 function
    FIG10 = task10()

    # Run task 11 function
    FIG11, focal_point = task11() #FOCAL_POINT = task11()

    # Run task 12 function
    FIG12 = task12()

    # Run task 13 function
    FIG13, TASK13_RMS = task13()

    # Run task 14 function
    # FIG14, TASK14_RMS, TASK14_DIFF_SCALE = task14()

    # Run task 15 function
    # FIG15_PC, FOCAL_POINT_PC, FIG15_CP, FOCAL_POINT_CP = task15()

    # Run task 16 function
    # FIG16, PC_RMS, CP_RMS, TASK16_DIFF_SCALE = task16()

    plt.show()
