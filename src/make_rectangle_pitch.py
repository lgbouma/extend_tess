import os
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import proj3d

def main():

    pointings = ['standard', 'c3p0', '54-40', '70-40', 'everyother']
    nsectors = [3, 4]
    elevs = [20, 90]

    for pointing in pointings:

        print(f"Starting {pointing}")

        for nsector in nsectors:

            for elev in elevs:

                # standard
                dl = 360/13
                if pointing == 'standard':
                    sc_elons = 120 + np.arange(dl - dl*nsector, dl, dl)[::-1]
                    sc_pitches = nsector * [54]
                    sc_erolls = nsector * [179.99]
                    view_kwargs = {'elev':elev, 'azim':90}

                # c3po
                elif pointing == 'c3p0':
                    sc_elons = 120 + np.arange(dl - dl*nsector, dl, dl)[::-1]
                    sc_pitches = nsector * [85]
                    sc_erolls = nsector * [179.99]
                    view_kwargs = {'elev':elev, 'azim':90}

                # 54/40
                elif pointing == '54-40':
                    elon0 = 40 if nsector == 3 else 50
                    sc_elons = elon0 + np.arange(2*dl - dl*nsector, 2*dl, dl)[::-1]
                    sc_pitches = nsector * [-54]
                    sc_erolls = nsector * [40]
                    view_kwargs = {'elev':elev, 'azim':90}

                elif pointing == '70-40':
                    elon0 = 40 if nsector == 3 else 50
                    sc_elons = elon0 + np.arange(2*dl - dl*nsector, 2*dl, dl)[::-1]
                    sc_pitches = nsector * [-70]
                    sc_erolls = nsector * [40]
                    view_kwargs = {'elev':elev, 'azim':90}

                elif pointing == 'everyother':
                    sc_elons = 195 + np.arange(dl - 2*dl*nsector, dl, 2*dl)[::-1]
                    sc_pitches = nsector * [54]
                    sc_erolls = nsector * [179.99]
                    view_kwargs = {'elev':elev, 'azim':90}

                plot_rectangles_on_sphere(pointing, sc_elons, sc_pitches, sc_erolls,
                                          nsector,
                                          view_kwargs=view_kwargs)



def plot_rectangles_on_sphere(
    pointing, sc_elons, sc_pitches, sc_erolls, nsector, view_kwargs={}
    ):

    fig = plt.figure(figsize=(4, 4), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # Plot the unit sphere
    u, v = np.mgrid[0:2*np.pi:500j, 0:np.pi:250j]
    r = 1.0
    x = r*np.cos(u)*np.sin(v)
    y = r*np.sin(u)*np.sin(v)
    z = r*np.cos(v)
    ax.plot_surface(x, y, z, color='white', alpha=0.05, linewidth=0)

    colors = [
        '#FFFF00',  # Yellow
        '#00FFFF',  # Cyan
        '#FF69B4',  # Hot Pink
        '#CCFF00',  # Electric Lime
        '#FF5E00',  # Bright Orange
        '#FF00FF',  # Magenta
        '#1F51FF',  # Neon Blue
    ]

    for idx, (elon, pitch, eroll) in enumerate(zip(sc_elons, sc_pitches, sc_erolls)):
        color = colors[idx % len(colors)]

        # Convert boresight coordinates to radians
        lon0 = np.radians(elon)
        pitch_rad = np.radians(pitch)
        roll = np.radians(eroll)

        # Initial boresight vector at latitude = 0
        n_boresight = np.array([
            np.cos(lon0),
            np.sin(lon0),
            0
        ])

        # Initial 'up' vector pointing north
        n_up = np.array([0, 0, 1])

        # Right vector orthogonal to boresight and up vectors
        n_right = np.cross(n_boresight, n_up)
        n_right = n_right / np.linalg.norm(n_right)

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

        # Rotate n_boresight and n_up_rot by pitch angle around n_right_rot
        R_pitch = rotation_matrix(n_right_rot, pitch_rad)
        n_boresight_rot = R_pitch.dot(n_boresight)
        n_up_rot2 = R_pitch.dot(n_up_rot)
        # n_right_rot remains the same

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
            delta = s * n_up_rot2 + t * n_right_rot
            delta_norm = np.linalg.norm(delta)
            if delta_norm != 0:
                axis = delta / delta_norm
                theta = delta_norm
                # Rotate n_boresight_rot by theta around axis
                R = rotation_matrix(axis, theta)
                n_point = R.dot(n_boresight_rot)
            else:
                n_point = n_boresight_rot
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

    # Plot equator
    max_theta = np.pi
    if 'elev' in view_kwargs:
        elev = view_kwargs['elev']
        if elev == 90:
            max_theta *= 2

    theta = np.linspace(0, max_theta, 100)
    x_eq = np.cos(theta)
    y_eq = np.sin(theta)
    z_eq = np.zeros_like(theta)
    ax.plot(x_eq, y_eq, z_eq, color='black', linewidth=0.5, ls=':')

    # Adjust the view
    ax.set_box_aspect([1, 1, 815/893])  # Equal aspect ratio
    ax.axis('off')

    # Apply view kwargs
    if 'elev' in view_kwargs and 'azim' in view_kwargs:
        ax.view_init(elev=view_kwargs['elev'], azim=view_kwargs['azim'])
        elev = view_kwargs['elev']
    else:
        ax.view_init(elev=20, azim=90)
        elev = 20

    outdir = '../results/EM3_rectangle_pitch_plot/'
    if not os.path.exists(outdir): os.mkdir(outdir)
    plt.savefig(join(outdir, f'{pointing}_N{nsector}_elev{elev}.pdf'))


if __name__ == "__main__":
    main()
