import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def generate_half_dome_waypoints(center, radius, steps, orientation):
    """
    Generate 3D waypoints forming a half-dome over a given center point, maintaining orientation.
    The points are distributed uniformly across the dome surface using an
    equal-area projection for the inclination angle.

    Parameters:
        center (tuple): (x, y, z) coordinates of the dome's base center.
        radius (float): Radius of the dome in meters.
        steps (int): Resolution of the dome (more steps = more waypoints).
        orientation (tuple): (rx, ry, rz) orientation to maintain for all waypoints.

    Returns:
        list of tuples: Waypoints as (x, y, z, rx, ry, rz) coordinates in meters and radians.
    """
    cx, cy, cz = center
    rx, ry, rz = orientation
    waypoints = []

    # Distribute points uniformly across the dome surface by spacing
    # phi using an equal area projection. This avoids clustering
    # of points near the top of the dome.
    for i in range(steps):
        # Generate phi so that the surface area band represented by
        # each step is approximately equal. The 0.5 offset prevents
        # points exactly at the pole or equator.
        phi = np.arccos(1 - (i + 0.5) / steps)

        # Theta is still sampled uniformly around the circumference.
        for j in range(steps * 2 + 1):  # More points around circumference
            theta = (2 * np.pi) * (j / (steps * 2))
            x = cx + radius * np.sin(phi) * np.cos(theta)
            y = cy + radius * np.sin(phi) * np.sin(theta)
            z = cz + radius * np.cos(phi)
            waypoints.append((x, y, z, rx, ry, rz))

    # Calculate the top of the dome position
    top_of_dome = (cx, cy, cz + radius)
    print(f"Center : {center}")
    print(f"Top of the Dome Position: {top_of_dome}")

    return waypoints

#Home Position: X: -0.170090, Y: -0.350146, Z: 0.108982, Rotation RX: -0.026828, RY: -0.005072, RZ: 1.497419 
#updated home position: X: -0.029437, Y: -0.295145, Z: 0.155930, Rotation RX: 0.039649, RY: 0.075470, RZ: -1.624816
# Update the function call to include orientation
points = generate_half_dome_waypoints(center=(-0.029437, -0.295145, 0.155930), radius=0.05, steps=15, orientation=(0.039674, 0.075407, -1.624840))

#python printout
print("[")
for each in points:
    print(f"    ({each[0]:.3f}, {each[1]:.3f}, {each[2]:.3f}, {each[3]:.3f}, {each[4]:.3f}, {each[5]:.3f}),")
print("]")

#cpp printout
for each in points:
    print(f"{{{each[0]:.3f}, {each[1]:.3f}, {each[2]:.3f}, {each[3]:.3f}, {each[4]:.3f}, {each[5]:.3f}}},")

# Extract x, y, z coordinates from points
x_coords, y_coords, z_coords = zip(*[(x, y, z) for x, y, z, _, _, _ in points])
center_x = sum(x_coords) / len(x_coords)
center_y = sum(y_coords) / len(y_coords)
center_z = sum(z_coords) / len(z_coords)

center_point = (center_x, center_y, center_z)
print(f"Calculated Center Point: {center_point}")

# Generate a color gradient based on the order of points
colors = np.linspace(0, 1, len(points))  # Create a gradient from 0 to 1

# Create a figure with 4 subplots
fig = plt.figure(figsize=(12, 10))

xmin = -0.3
xmax = 0.3
ymin = -0.3
ymax = 0.3
zmin = 0.0
zmax = 0.3

# 3D view
ax1 = fig.add_subplot(221, projection='3d')
scatter = ax1.scatter(x_coords, y_coords, z_coords, c=colors, cmap='viridis', marker='o')
ax1.set_xlim([xmin, xmax])
ax1.set_ylim([ymin, ymax])
ax1.set_zlim([zmin, zmax])
ax1.set_title("3D View (Point Order)")
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')
ax1.grid(True)
fig.colorbar(scatter, ax=ax1, label="Point Order")

# Top view (XY plane)
ax2 = fig.add_subplot(222)
scatter = ax2.scatter(x_coords, y_coords, c=colors, cmap='viridis', marker='o')
ax2.set_xlim([xmin, xmax])
ax2.set_ylim([ymin, ymax])
ax2.set_title("Top View (XY Plane, Point Order)")
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.grid(True)
fig.colorbar(scatter, ax=ax2, label="Point Order")

# Front view (XZ plane)
ax3 = fig.add_subplot(223)
scatter = ax3.scatter(x_coords, z_coords, c=colors, cmap='viridis', marker='o')
ax3.set_xlim([xmin, xmax])
ax3.set_ylim([zmin, zmax])
ax3.set_title("Front View (XZ Plane, Point Order)")
ax3.set_xlabel('X')
ax3.set_ylabel('Z')
ax3.grid(True)
fig.colorbar(scatter, ax=ax3, label="Point Order")

# Side view (YZ plane)
ax4 = fig.add_subplot(224)
scatter = ax4.scatter(y_coords, z_coords, c=colors, cmap='viridis', marker='o')
ax4.set_xlim([ymin, ymax])
ax4.set_ylim([zmin, zmax])
ax4.set_title("Side View (YZ Plane, Point Order)")
ax4.set_xlabel('Y')
ax4.set_ylabel('Z')
ax4.grid(True)
fig.colorbar(scatter, ax=ax4, label="Point Order")

# Adjust layout and show the plot
plt.tight_layout()
plt.show()
