import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import proj3d

def plot_rectangles_on_sphere(pointing, sc_elons, sc_elats, sc_erolls, view_kwargs={}):
    fig = plt.figure(figsize=(4, 4), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # Plot the unit sphere
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    ax.plot_surface(x, y, z, color='w', alpha=0.05, linewidth=0)

    colors = [
        '#FFFF00',  # Yellow
        '#00FFFF',  # Cyan
        '#FF69B4',  # Hot Pink
        '#CCFF00',  # Electric Lime
        '#FF5E00',  # Bright Orange
        '#FF00FF',  # Magenta
        '#1F51FF',  # Neon Blue
    ]

    for idx, (elon, elat, eroll) in enumerate(zip(sc_elons, sc_elats, sc_erolls)):
        color = colors[idx % len(colors)]

        # Convert boresight coordinates to radians
        lon0 = np.radians(elon)
        lat0 = np.radians(elat)
        roll = np.radians(eroll)

        # Compute boresight vector
        n_boresight = np.array([
            np.cos(lat0) * np.cos(lon0),
            np.cos(lat0) * np.sin(lon0),
            np.sin(lat0)
        ])

        # Compute local coordinate axes
        n_up = np.array([
            -np.sin(lat0) * np.cos(lon0),
            -np.sin(lat0) * np.sin(lon0),
            np.cos(lat0)
        ])

        n_right = np.array([
            -np.sin(lon0),
            np.cos(lon0),
            0
        ])

        # Rotate n_up and n_right by roll angle around n_boresight
        def rotation_matrix(axis, theta):
            axis = axis / np.linalg.norm(axis)
            a = np.cos(theta / 2.0)
            b, c, d = -axis * np.sin(theta / 2.0)
            return np.array([
                [a*a + b*b - c*c - d*d, 2*(b*c - a*d),     2*(b*d + a*c)],
                [2*(b*c + a*d),     a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
                [2*(b*d - a*c),     2*(c*d + a*b),     a*a + d*d - b*b - c*c]
            ])

        R_roll = rotation_matrix(n_boresight, roll)
        n_up_rot = R_roll.dot(n_up)
        n_right_rot = R_roll.dot(n_right)

        # Create grid over the rectangle
        L_up = np.radians(12)  # Half of 24 degrees
        L_right = np.radians(48)  # Half of 96 degrees
        s_vals = np.linspace(-L_up, L_up, 10)
        t_vals = np.linspace(-L_right, L_right, 20)
        s_grid, t_grid = np.meshgrid(s_vals, t_vals)
        s_flat = s_grid.flatten()
        t_flat = t_grid.flatten()

        # For each point, compute rotation axis and angle
        n_points = []
        for s, t in zip(s_flat, t_flat):
            delta = s * n_up_rot + t * n_right_rot
            delta_norm = np.linalg.norm(delta)
            if delta_norm != 0:
                axis = delta / delta_norm
                theta = delta_norm
                # Rotate n_boresight by theta around axis
                R = rotation_matrix(axis, theta)
                n_point = R.dot(n_boresight)
            else:
                n_point = n_boresight
            n_points.append(n_point)

        n_points = np.array(n_points)

        x_points = n_points[:, 0]
        y_points = n_points[:, 1]
        z_points = n_points[:, 2]

        # Reshape back to grid shape
        x_grid = x_points.reshape(s_grid.shape)
        y_grid = y_points.reshape(s_grid.shape)
        z_grid = z_points.reshape(s_grid.shape)

        # Plot the surface
        ax.plot_surface(
            x_grid, y_grid, z_grid,
            color=color, alpha=0.3,
            edgecolor='black', linewidth=0.25,
            rcount=4, ccount=1
        )

        # Plot the border
        # Extract the perimeter points
        border_points = []
        # Top edge
        border_points.extend(zip(x_grid[0, :], y_grid[0, :], z_grid[0, :]))
        # Right edge
        border_points.extend(zip(x_grid[:, -1], y_grid[:, -1], z_grid[:, -1]))
        # Bottom edge (reversed)
        border_points.extend(zip(x_grid[-1, ::-1], y_grid[-1, ::-1], z_grid[-1, ::-1]))
        # Left edge (reversed)
        border_points.extend(zip(x_grid[::-1, 0], y_grid[::-1, 0], z_grid[::-1, 0]))
        border_points = np.array(border_points)

        ax.plot(
            border_points[:, 0], border_points[:, 1], border_points[:, 2],
            color='black', linewidth=0.5
        )

    # Add the arrow from the center to slightly beyond the north pole
    ax.quiver(
        0, 0, 1,  # Starting point
        0, 0, 1.2,  # Direction vector
        length=0.2, color='black', alpha=0.5,
        arrow_length_ratio=0.4, linewidth=0.5
    )

    # plot equator
    theta = np.linspace(0, 2 * np.pi, 100)
    x_eq = np.cos(theta)
    y_eq = np.sin(theta)
    z_eq = np.zeros_like(theta)
    ax.plot(x_eq, y_eq, z_eq, color='black', linewidth=0.5, ls=':')

    # Adjust the view
    ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
    ax.axis('off')

    # Apply view kwargs
    if 'elev' in view_kwargs and 'azim' in view_kwargs:
        ax.view_init(elev=view_kwargs['elev'], azim=view_kwargs['azim'])
    else:
        ax.view_init(elev=20, azim=90)

    plt.savefig(f'../results/EM3_rectangle_sphere_plot/{pointing}.pdf')

# Example usage

pointings = ['standard', 'c3p0', '54-40']

for pointing in pointings:

    # standard
    dl = 360/13
    if pointing == 'standard':
        sc_elons = 120 + np.array([dl, 0, -dl])
        sc_elats = [54, 54, 54]
        sc_erolls = [179.99, 179.99, 179.99]
        view_kwargs = {'elev':20, 'azim':90}

    # c3po
    elif pointing == 'c3p0':
        sc_elons = 120 + np.array([dl, 0, -dl])
        sc_elats = [85, 85, 85]
        sc_erolls = [179.99, 179.99, 179.99]
        view_kwargs = {'elev':20, 'azim':90}

    # 54/40
    elif pointing == '54-40':
        sc_elons = 100 + np.array([dl, 0, -dl, -2*dl])
        sc_elats = [54, 54, 54, 54]
        sc_erolls = [40, 40, 40, 40]
        view_kwargs = {'elev':20, 'azim':90}

    plot_rectangles_on_sphere(pointing, sc_elons, sc_elats, sc_erolls,
                              view_kwargs=view_kwargs)
