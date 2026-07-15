# import library 
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import differential_evolution
import pandas as pd
from plyfile import PlyData, PlyElement

# =============================================================================
# (Unchanged) Functions from your original code:
# =============================================================================

def xi(x):
    return 0.5 * (1 - x / np.sqrt(1 + x * x))


def beta_equ(x, P, Rs):
    return x ** P[0] * Rs[0] + x ** P[1] * Rs[1] + x ** P[2] * Rs[2] - 1


def calculate_triangle_centers(vertices, faces):
    """
    Calculate the incenter of each triangle in a triangulation.
    """
    triangle_centers = []  # Initialize an empty list to store triangle centers
    for face in faces:  # Iterate through each triangle's vertex indices
        v1 = vertices[face[0]]
        v2 = vertices[face[1]]
        v3 = vertices[face[2]]
        a = np.linalg.norm(v2 - v3)
        b = np.linalg.norm(v1 - v3)
        c = np.linalg.norm(v1 - v2)
        center = (a * v1 + b * v2 + c * v3) / (a + b + c)
        triangle_centers.append(center)
    return np.array(triangle_centers)


def calculate_triangle_areas(vertices, faces):
    """
    Calculate the area of each triangle in a triangulation.
    """
    triangle_areas = []
    for face in faces:
        v1 = vertices[face[0]]
        v2 = vertices[face[1]]
        v3 = vertices[face[2]]
        a = v2 - v1
        b = v3 - v1
        area = 0.5 * np.linalg.norm(np.cross(a, b))
        triangle_areas.append(area)
    return np.array(triangle_areas)


def calculate_triangle_normals(center_points, cell_radii, cell_shape_factor):
    """
    Calculate the normals to the center of each triangle.
    """
    nx = cell_shape_factor[0] * np.sign(center_points[:, 0]) * \
         np.abs(center_points[:, 0]) ** (cell_shape_factor[0] - 1) / (cell_radii[0] ** cell_shape_factor[0])
    ny = cell_shape_factor[1] * np.sign(center_points[:, 1]) * \
         np.abs(center_points[:, 1]) ** (cell_shape_factor[1] - 1) / (cell_radii[1] ** cell_shape_factor[1])
    nz = cell_shape_factor[2] * np.sign(center_points[:, 2]) * \
         np.abs(center_points[:, 2]) ** (cell_shape_factor[2] - 1) / (cell_radii[2] ** cell_shape_factor[2])
    normals = np.column_stack((nx, ny, nz))
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)
    return normals


def project_to_superellipsoid(center_points, cell_radii, cell_shape_factor):
    """
    Project triangle centers to the surface of a superellipsoid.
    """
    projected_centers = np.zeros_like(center_points)
    for i in range(len(center_points)):
        c = center_points[i]
        cs = (np.abs(c) / cell_radii) ** cell_shape_factor
        a, b = 0, 1
        while beta_equ(a, cell_shape_factor, cs) * beta_equ(b, cell_shape_factor, cs) > 0:
            b *= 2
        xs = root_scalar(lambda x: beta_equ(x, cell_shape_factor, cs), bracket=[a, b]).root
        projected_centers[i] = xs * c
    return projected_centers


def load_triangulated_surface(file_path):
    """
    Load a triangulated surface from a .mat file.
    """
    mat_data = scipy.io.loadmat(file_path)
    vertices = mat_data['V']
    faces = mat_data['F'] - 1  # Adjust face indices to start from 0
    return vertices, faces


def process_triangulated_surface(file_path, cell_radii, cell_shape_factor):
    """
    Process a triangulated surface and return projected centers, normals, triangle areas, etc.
    """
    vertices, faces = load_triangulated_surface(file_path)
    vertices = vertices * cell_radii
    center_points = calculate_triangle_centers(vertices, faces)
    projected_centers = project_to_superellipsoid(center_points, cell_radii, cell_shape_factor)
    normals = calculate_triangle_normals(projected_centers, cell_radii, cell_shape_factor)
    triangle_areas = calculate_triangle_areas(vertices, faces)
    num_triangles = len(center_points)
    total_area = np.sum(triangle_areas)
    return projected_centers, normals, triangle_areas, num_triangles, total_area


def save_cell_bindings_to_file(centrosomes, filename):
    """
    Save cell bindings to a file in multiple columns.

    Parameters:
        centrosomes (dict): Dictionary containing the binding probabilities for centrosomes with the cell.
        filename (str): Path to the output file.
    """
    with open(filename, 'w') as f:
        f.write("cent1_cell_binding,cent2_cell_binding\n")
        cell_binding = np.column_stack((centrosomes["cent1_cell_binding"], centrosomes["cent2_cell_binding"]))
        np.savetxt(f, cell_binding, delimiter=',', fmt='%0.8f')

def save_surface_ply(output_ply, vertices, faces, cent_p, min_val, max_val):
    """
    Save the full triangulated surface as a PLY file with color-coded binding values.

    Parameters:
        output_ply (str): Path to save the output `.ply` file.
        vertices (ndarray): (N, 3) array of vertex positions.
        faces (ndarray): (M, 3) array of triangle indices.
        cent_p (ndarray): (M, 2) array of binding values for each triangle.
        min_val (float): Lower threshold for binding values (black if below this).
        max_val (float): Upper threshold for binding values (white if above this).
    """

    # Compute the average binding values per triangle
    binding_value = (cent_p[:, 0] + cent_p[:, 1]) / 2

    # Interpolate binding values to original vertices
    vertex_binding = np.zeros(vertices.shape[0])  # Initialize per-vertex binding values
    vertex_counts = np.zeros(vertices.shape[0])   # Count how many triangles each vertex belongs to

    # Assign each vertex the average binding value of the triangles it belongs to
    for i, face in enumerate(faces):
        for v in face:  # Assign binding value to each vertex of the triangle
            vertex_binding[v] += binding_value[i]
            vertex_counts[v] += 1

    # Normalize binding values by averaging
    vertex_binding /= np.maximum(vertex_counts, 1)  # Avoid division by zero

    # Normalize for grayscale mapping
    vertex_colors = np.zeros_like(vertex_binding)  # Default: black
    mask = (vertex_binding >= min_val) & (vertex_binding <= max_val)

    if max_val > min_val:
        vertex_colors[mask] = 255 * (vertex_binding[mask] - min_val) / (max_val - min_val)

    vertex_colors[vertex_binding > max_val] = 255  # White if above max_val

    # Convert to uint8 for PLY vertex colors
    colors = np.column_stack((vertex_colors, vertex_colors, vertex_colors)).astype(np.uint8)

    # Prepare vertex elements
    vertex_data = np.array([
        (vertices[i, 0], vertices[i, 1], vertices[i, 2], 
         colors[i, 0], colors[i, 1], colors[i, 2]) for i in range(len(vertices))
    ], dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('red', 'u1'), ('green', 'u1'), ('blue', 'u1')])

    vertex_element = PlyElement.describe(vertex_data, 'vertex')

    # Prepare face elements
    face_data = np.array([(list(face),) for face in faces], dtype=[('vertex_indices', 'i4', (3,))])
    face_element = PlyElement.describe(face_data, 'face')

    # Save PLY file
    ply_data = PlyData([vertex_element, face_element], text=True)
    ply_data.write(output_ply)
    print(f"PLY file saved as {output_ply} (Surface with Colors)")

def visualize_surface(vertices, faces):
    """
    Visualize a triangulated surface using matplotlib.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(vertices[:, 0], vertices[:, 1], vertices[:, 2], triangles=faces)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Triangulated Surface')
    max_range = np.array([vertices[:, 0].max() - vertices[:, 0].min(),
                          vertices[:, 1].max() - vertices[:, 1].min(),
                          vertices[:, 2].max() - vertices[:, 2].min()]).max() / 2.0
    mean_x = vertices[:, 0].mean()
    mean_y = vertices[:, 1].mean()
    mean_z = vertices[:, 2].mean()
    ax.set_xlim(mean_x - max_range, mean_x + max_range)
    ax.set_ylim(mean_y - max_range, mean_y + max_range)
    ax.set_zlim(mean_z - max_range, mean_z + max_range)
    ax.set_box_aspect([1, 1, 1])
    plt.show()


def update_kesi(cent_pos, triangle_center):
    """
    Update the normalized vectors (kesi) and distances (D) from centrosomes to surface elements.
    """
    cent1_kesi = triangle_center - cent_pos[0]
    cent2_kesi = triangle_center - cent_pos[1]
    cent1_d = np.linalg.norm(cent1_kesi, axis=1)
    cent2_d = np.linalg.norm(cent2_kesi, axis=1)
    cent1_kesi /= cent1_d[:, np.newaxis]
    cent2_kesi /= cent2_d[:, np.newaxis]
    cent_kesi = np.hstack((cent1_kesi, cent2_kesi))
    cent_d = np.hstack((cent1_d[:, np.newaxis], cent2_d[:, np.newaxis]))
    return cent_kesi, cent_d


def update_force(cent_kesi, cent_p, triangle_area, total_area, cfg_m, cfg_f0):
    """
    Update the forces on centrosomes from cortical force-generators.
    """
    alpha = cfg_m * cfg_f0 / total_area
    cent1_f = np.column_stack((
        alpha * cent_p[:, 0] * cent_kesi[:, 0] * triangle_area,
        alpha * cent_p[:, 0] * cent_kesi[:, 1] * triangle_area,
        alpha * cent_p[:, 0] * cent_kesi[:, 2] * triangle_area
    ))
    cent2_f = np.column_stack((
        alpha * cent_p[:, 1] * cent_kesi[:, 3] * triangle_area,
        alpha * cent_p[:, 1] * cent_kesi[:, 4] * triangle_area,
        alpha * cent_p[:, 1] * cent_kesi[:, 5] * triangle_area
    ))
    cent_f = np.array([np.sum(cent1_f, axis=0), np.sum(cent2_f, axis=0)])
    return cent_f


def update_cent_pos(cent_f, cent_pos, eta, nu, l_dot, deltat):
    """
    Update the position and velocity of centrosomes.
    """
    r1 = cent_pos[0, :]
    r2 = cent_pos[1, :]
    l = np.linalg.norm(r1 - r2)
    s = (r1 - r2)
    s_hat = s / l
    f1 = cent_f[0, :]
    f2 = cent_f[1, :]
    df = f1 - f2 + nu * l_dot * s_hat
    ft = f1 + f2

    vc = ft / (2.0 * eta)
    vs = 0.5 * np.dot(df, s_hat) * s_hat / (eta + nu) + 0.5 * df / eta - 0.5 * np.dot(df, s_hat) * s_hat / eta

    v1 = vc + vs
    v2 = vc - vs

    r1 = r1 + v1 * deltat
    r2 = r2 + v2 * deltat

    cent_vel = np.array([v1, v2])
    cent_pos = np.array([r1, r2])
    return cent_pos, cent_vel


def update_impingement(mt_vg, mt_nuc_rate, mt_ave_length, cfg_r, triangle_center, triangle_normal, cent_vel, cent_pos):
    """
    Update the impingement rate of microtubules to surface elements.
    """
    alpha = mt_nuc_rate / mt_vg
    # For centrosome 1
    x1 = triangle_center[:, 0] - cent_pos[0, 0]
    x2 = triangle_center[:, 1] - cent_pos[0, 1]
    x3 = triangle_center[:, 2] - cent_pos[0, 2]
    d1 = np.sqrt(x1**2 + x2**2 + x3**2)
    kesi1x = x1 / d1
    kesi1y = x2 / d1
    kesi1z = x3 / d1
    impinge1 = (cent_vel[0, 0] + mt_vg * kesi1x) * triangle_normal[:, 0] + \
               (cent_vel[0, 1] + mt_vg * kesi1y) * triangle_normal[:, 1] + \
               (cent_vel[0, 2] + mt_vg * kesi1z) * triangle_normal[:, 2]
    impinge1[impinge1 < 0] = 0
    omega1 = alpha * xi(d1 / cfg_r) * np.exp(-d1 / mt_ave_length) * impinge1
    # For centrosome 2
    x1 = triangle_center[:, 0] - cent_pos[1, 0]
    x2 = triangle_center[:, 1] - cent_pos[1, 1]
    x3 = triangle_center[:, 2] - cent_pos[1, 2]
    d2 = np.sqrt(x1**2 + x2**2 + x3**2)
    kesi2x = x1 / d2
    kesi2y = x2 / d2
    kesi2z = x3 / d2
    impinge2 = (cent_vel[1, 0] + mt_vg * kesi2x) * triangle_normal[:, 0] + \
               (cent_vel[1, 1] + mt_vg * kesi2y) * triangle_normal[:, 1] + \
               (cent_vel[1, 2] + mt_vg * kesi2z) * triangle_normal[:, 2]
    impinge2[impinge2 < 0] = 0
    omega2 = alpha * xi(d2 / cfg_r) * np.exp(-d2 / mt_ave_length) * impinge2
    return omega1, omega2


def update_p(mt_omega1, mt_omega2, cent_p, cfg_m, cfg_m0, cfg_k, deltat):
    """
    Update the probability of attachment of microtubules to cortical force-generators.
    """
    alpha = cfg_m / cfg_m0
    alpha_inv = 1 / alpha
    om1 = mt_omega1 * alpha_inv
    om2 = mt_omega2 * alpha_inv
    exp_term = np.exp(-alpha * (1 - cent_p[:, 0] - cent_p[:, 1]))
    delta_p1 = om1 * (1 - exp_term) - cfg_k * cent_p[:, 0]
    delta_p2 = om2 * (1 - exp_term) - cfg_k * cent_p[:, 1]
    cent_p += np.column_stack((delta_p1, delta_p2)) * deltat
    cent_p[cent_p[:, 0] < 0, 0] = 0
    cent_p[cent_p[:, 1] < 0, 1] = 0
    return cent_p

# =============================================================================
# New Functions for Optimization:
# =============================================================================

def simulate_simulation(params):
    """
    Run the simulation with a given set of parameters.
    The parameters to optimize are:
      params[0]: mt_inf
      params[1]: mt_vg
      params[2]: mt_ave_length
      params[3]: cfg_k
      params[4]: cfg_m
      params[5]: eta
      params[6]: nu
    Returns the simulation history (cent_pos_array), where:
      - Column 0 is time,
      - Column 7 is the spindle length (distance between centrosomes).
    """
    # Unpack optimized parameters:
    mt_inf, mt_vg, mt_ave_length, cfg_k, cfg_m, eta, nu = params
    # Fixed parameters:
    cfg_f0 = 10
    cfg_r = 1.0
    deltat = 0.1
    t_max = 200.0
    initial_spindle_length = 13.26534077
    cell_length = 22.5
    mt_cat_rate = mt_vg / mt_ave_length
    mt_nuc_rate = mt_inf * mt_cat_rate
    chrom_a = 6.23
    chrom_T = 28.7

    # Load and process the cell surface:
    file_path = './unit_sphere.mat'
    cell_radii = np.array([cell_length, cell_length, cell_length])
    cell_shape_factor = np.array([2.0, 2.0, 2.0])
    triangle_center, triangle_normal, triangle_area, num_triangles, total_area = process_triangulated_surface(
        file_path, cell_radii, cell_shape_factor)

    cfg_m0 = total_area / (np.pi * cfg_r * cfg_r)
    cent_pos = np.array([[initial_spindle_length / 2, 0, 0],
                         [-initial_spindle_length / 2, 0, 0]])
    cent_p = np.zeros((num_triangles, 2))
    cent_pos_array = []
    cent_vel = np.array([[0, 0, 0], [0, 0, 0]])

    # Pre-equilibrate microtubule attachment probabilities:
    omega1, omega2 = update_impingement(mt_vg, mt_nuc_rate, mt_ave_length, cfg_r,
                                        triangle_center, triangle_normal, cent_vel, cent_pos)
    for _ in range(10000):
        cent_p = update_p(omega1, omega2, cent_p, cfg_m, cfg_m0, cfg_k, deltat)

    # Simulation time-stepping loop:
    for t in np.arange(0, t_max + deltat, deltat):
        l_dot = chrom_a * np.exp(-t / chrom_T) / chrom_T
        cent_kesi, cent_d = update_kesi(cent_pos, triangle_center)
        cent_f = update_force(cent_kesi, cent_p, triangle_area, total_area, cfg_m, cfg_f0)
        cent_pos, cent_vel = update_cent_pos(cent_f, cent_pos, eta, nu, l_dot, deltat)
        omega1, omega2 = update_impingement(mt_vg, mt_nuc_rate, mt_ave_length, cfg_r,
                                            triangle_center, triangle_normal, cent_vel, cent_pos)
        cent_p = update_p(omega1, omega2, cent_p, cfg_m, cfg_m0, cfg_k, deltat)
        cent_dist = np.linalg.norm(cent_pos[0, :] - cent_pos[1, :])
        # Record: time, centrosome positions, and spindle length (cent_dist)
        cent_pos_array.append(np.concatenate(([t], cent_pos[0, :], cent_pos[1, :], [cent_dist])))
    cent_pos_array = np.array(cent_pos_array)
    return cent_pos_array


def simulation_error(params, exp_time, exp_spindle):
    """
    Given a parameter vector and experimental data (exp_time and exp_spindle),
    run the simulation and compute the error as the sum of squared differences
    between the simulated spindle length (column 7) and the experimental spindle length.
    """
    sim_data = simulate_simulation(params)
    sim_time = sim_data[:, 0]
    sim_spindle = sim_data[:, 7]
    # Interpolate simulation spindle length at experimental time points:
    sim_spindle_interp = np.interp(exp_time, sim_time, sim_spindle)
    error = np.sum((sim_spindle_interp - exp_spindle) ** 2)
    return error

# =============================================================================
# Main Function: Optimization and Plotting
# =============================================================================

def main_custom():
    # Define simulation parameters
    mt_inf = 5000
    cfg_f0 = 5
    cfg_r = 1.5
    cfg_k = 0.05
    mt_vg = 1
    eta = 450
    nu = 50
    mt_cat_rate = 0.03839612
    deltat = 0.1
    t_max = 206.0
    cfg_m = 200
    initial_spindle_length = 13
    cell_length = 23
    chrom_a = 6.23
    chrom_T = 28.7
    mt_ave_length = mt_vg / mt_cat_rate
    mt_nuc_rate = mt_inf * mt_cat_rate

    # Load the cell surface
    file_path = './unit_sphere.mat'
    cell_radii = np.array([cell_length, cell_length, cell_length])
    vertices, faces = load_triangulated_surface(file_path)
    cell_shape_factor = np.array([2.0, 2.0, 2.0])
    triangle_center, triangle_normal, triangle_area, num_triangles, total_area = \
        process_triangulated_surface(file_path, cell_radii, cell_shape_factor)

    # Initialize simulation state
    cfg_m0 = total_area / (np.pi * cfg_r * cfg_r)
    cent_pos = np.array([[initial_spindle_length / 2, 0, 0],
                         [-initial_spindle_length / 2, 0, 0]])
    cent_p = np.ones((num_triangles, 2))
    cent_pos_array = []
    cent_vel = np.array([[0, 0, 0], [0, 0, 0]])

    # Pre-equilibrate microtubule attachment probabilities
    omega1, omega2 = update_impingement(mt_vg, mt_nuc_rate, mt_ave_length, cfg_r, 
                                        triangle_center, triangle_normal, cent_vel, cent_pos)
    for _ in range(10000):
        cent_p = update_p(omega1, omega2, cent_p, cfg_m, cfg_m0, cfg_k, deltat)

    # Run simulation loop
    for t in np.arange(0, t_max + deltat, deltat):
        l_dot = chrom_a * np.exp(-t/chrom_T) / chrom_T
        cent_kesi, cent_d = update_kesi(cent_pos, triangle_center)
        cent_f = update_force(cent_kesi, cent_p, triangle_area, total_area, cfg_m, cfg_f0)
        cent_pos, cent_vel = update_cent_pos(cent_f, cent_pos, eta, nu, l_dot, deltat)

        # Update impingement rate and probabilities
        omega1, omega2 = update_impingement(mt_vg, mt_nuc_rate, mt_ave_length, cfg_r, 
                                            triangle_center, triangle_normal, cent_vel, cent_pos)
        cent_p = update_p(omega1, omega2, cent_p, cfg_m, cfg_m0, cfg_k, deltat)

        # Store centrosome positions and spindle length
        cent_dist = np.linalg.norm(cent_pos[0, :] - cent_pos[1, :])
        cent_pos_array.append(np.concatenate(([t], cent_pos[0, :], cent_pos[1, :], [cent_dist])))
        

    cent_pos_array = np.array(cent_pos_array)

    # Load experimental data
    exp_data = np.loadtxt('2cell.txt')
    exp_time = exp_data[:, 0]
    exp_spindle = exp_data[:, 1]

    # Plot results
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Plot centrosome X positions vs. Time
    axes[0].plot(cent_pos_array[:, 0], cent_pos_array[:, 1], label='Centrosome 1')
    axes[0].plot(cent_pos_array[:, 0], cent_pos_array[:, 4], label='Centrosome 2')
    axes[0].set_xlabel('Time [s]')
    axes[0].set_ylabel('Centrosome X [micron]')
    axes[0].set_title('Centrosome X Position vs. Time')
    axes[0].grid(True)
    axes[0].legend()

    # Plot spindle length vs. Time (Simulation vs Experiment)
    axes[1].plot(cent_pos_array[:, 0], cent_pos_array[:, 7], label='Simulation', color='blue')
    axes[1].scatter(exp_time, exp_spindle, label='Experiment', color='red', marker='o', s=40)
    axes[1].set_xlabel('Time [s]')
    axes[1].set_ylabel('Spindle Length [micron]')
    axes[1].set_title('Spindle Length vs. Time')
    axes[1].grid(True)
    axes[1].legend()
    # Create a DataFrame
    df = pd.DataFrame({
        'Time [s]': cent_pos_array[:, 0],
        'Spindle Length [micron]': cent_pos_array[:, 7]
    })

    plt.tight_layout()
    plt.show()

    plt.tight_layout()
    output_filename = "simulation_results_2cell.png"
    df.to_csv('spindle_length_2cell.csv', index=False)
    plt.savefig(output_filename, dpi=300)
    plt.close(fig)
    print(f"Plot saved as {output_filename}")        

    save_surface_ply("./binding_2cell_t200.ply", vertices, faces, cent_p, 0.3, 0.7)

if __name__ == "__main__":
    main_custom()
