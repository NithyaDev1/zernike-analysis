import csv
import numpy as np
import matplotlib.pyplot as plt


# ==========================================================
# SETTINGS
# ==========================================================

folder = "20261013"

positions = [
    "P0", "P1", "P2", "P3", "P4",
    "P5", "P6", "P7", "P8", "P9"
]

# Source-to-model-cornea distance at P0
source_distance_P0 = 23.43   # mm

# Each backwards translation
step_size = 0.50             # mm


# ==========================================================
# OPTICAL SYSTEM VALUES
# ==========================================================

# Model cornea lens
f_cornea = 25.0

# Model cornea -> L1
d_cornea_L1 = 58.43

# L1
f1 = 46.2

# L1 -> L2
d_L1_L2 = 201.55

# L2
f2 = 75.0

# L2 -> L3
d_L2_L3 = 96.75

# L3
f3 = 100.0

# UPDATED L3 -> Shack-Hartmann sensor
d_L3_SH = 91.2


# ==========================================================
# SHACK-HARTMANN PUPIL
# ==========================================================

# Thorlabs pupil diameter = 3 mm
pupil_radius = 1.5   # mm


# ==========================================================
# ABCD MATRICES
# ==========================================================

def P(d):
    """
    Free-space propagation over distance d.
    """
    return np.array([
        [1.0, d],
        [0.0, 1.0]
    ])


def L(f):
    """
    Thin lens with focal length f.
    """
    return np.array([
        [1.0, 0.0],
        [-1.0 / f, 1.0]
    ])


# ==========================================================
# STORAGE
# ==========================================================

predicted_z5_raw = []
thorlabs_z5 = []
curvatures = []
matrices = []
source_distances = []


# ==========================================================
# LOOP THROUGH P0-P9
# ==========================================================

for i, position in enumerate(positions):

    # ------------------------------------------------------
    # Source distance for this position
    # ------------------------------------------------------

    d_source = source_distance_P0 + step_size * i

    source_distances.append(d_source)


    # ------------------------------------------------------
    # Build full ABCD matrix
    #
    # Source
    # -> cornea lens
    # -> L1
    # -> L2
    # -> L3
    # -> SH sensor
    # ------------------------------------------------------

    M = (
        P(d_L3_SH)
        @ L(f3)
        @ P(d_L2_L3)
        @ L(f2)
        @ P(d_L1_L2)
        @ L(f1)
        @ P(d_cornea_L1)
        @ L(f_cornea)
        @ P(d_source)
    )

    matrices.append(M)


    # ------------------------------------------------------
    # Extract ABCD values
    # ------------------------------------------------------

    A = M[0, 0]
    B = M[0, 1]
    C = M[1, 0]
    D = M[1, 1]


    # ------------------------------------------------------
    # Radius of curvature at SH plane
    #
    # For an on-axis point source:
    #
    # x_SH     = B * theta_0
    # theta_SH = D * theta_0
    #
    # Therefore:
    #
    # theta_SH / x_SH = D / B = 1/R
    #
    # so:
    #
    # R = B / D
    # ------------------------------------------------------

    R = B / D

    curvatures.append(R)


    # ------------------------------------------------------
    # Convert wavefront curvature into Z5
    #
    # W(r) = r^2 / (2R)
    #
    # Thorlabs Z5:
    #
    # Z5 = sqrt(3) * (2 rho^2 - 1)
    #
    # Therefore:
    #
    # C5 = a^2 / (4 sqrt(3) R)
    #
    # x1000 converts mm -> micrometers
    # ------------------------------------------------------

    z5 = (
        1000
        * pupil_radius**2
        /
        (
            4
            * np.sqrt(3)
            * R
        )
    )

    predicted_z5_raw.append(z5)


    # ------------------------------------------------------
    # Read Thorlabs Z5 from MD file
    #
    # Z5 is Excel row 49
    # Python row index = 48
    # coefficient is column 4 -> index 3
    # ------------------------------------------------------

    filename = f"{folder}/1013{position}MD.csv"

    with open(filename, "r", encoding="cp1252") as file:
        rows = list(csv.reader(file))

    z5_thorlabs = float(
        rows[48][3]
    )

    thorlabs_z5.append(
        z5_thorlabs
    )


# ==========================================================
# CONVERT TO NUMPY ARRAYS
# ==========================================================

predicted_z5_raw = np.array(
    predicted_z5_raw
)

thorlabs_z5 = np.array(
    thorlabs_z5
)

curvatures = np.array(
    curvatures
)


# ==========================================================
# COMPARE CHANGE RELATIVE TO P0
# ==========================================================

# Our ABCD calculation predicts that the raw curvature
# coefficient decreases as the source moves backward.
#
# Thorlabs uses the opposite sign direction in our
# experimental trend, so we flip the predicted CHANGE.
#
# This does NOT change the magnitude.
# ==========================================================

predicted_delta_z5 = -(
    predicted_z5_raw
    - predicted_z5_raw[0]
)

thorlabs_delta_z5 = (
    thorlabs_z5
    - thorlabs_z5[0]
)


# ==========================================================
# PRINT RESULTS
# ==========================================================

print()
print("PREDICTED Z5 vs THORLABS Z5")
print("=" * 83)

print(
    f"{'Pos':<6}"
    f"{'Source d':<12}"
    f"{'R at SH':<12}"
    f"{'Pred ΔZ5':<14}"
    f"{'Thor Z5':<14}"
    f"{'Thor ΔZ5':<14}"
)

print("-" * 83)


for i in range(len(positions)):

    print(
        f"{positions[i]:<6}"
        f"{source_distances[i]:<12.2f}"
        f"{curvatures[i]:<12.3f}"
        f"{predicted_delta_z5[i]:<14.4f}"
        f"{thorlabs_z5[i]:<14.4f}"
        f"{thorlabs_delta_z5[i]:<14.4f}"
    )


# ==========================================================
# PRINT EACH ABCD MATRIX
# ==========================================================

print()
print("ABCD MATRICES")
print("=" * 60)


for i, position in enumerate(positions):

    print()
    print(
        f"{position}   "
        f"source distance = "
        f"{source_distances[i]:.2f} mm"
    )

    print(
        matrices[i]
    )


# ==========================================================
# PLOT
# ==========================================================

x = np.arange(
    len(positions)
)


plt.figure(
    figsize=(9, 5)
)

plt.plot(
    x,
    predicted_delta_z5,
    marker="o",
    label="Predicted Z5 change (ABCD)"
)

plt.plot(
    x,
    thorlabs_delta_z5,
    marker="o",
    label="Measured Z5 change (Thorlabs)"
)

plt.xticks(
    x,
    positions
)

plt.xlabel(
    "Pinhole Position"
)

plt.ylabel(
    "Change in Z5 from P0 (µm)"
)

plt.title(
    "Predicted vs Measured Defocus"
)

plt.legend()

plt.grid(
    alpha=0.3
)

plt.tight_layout()

plt.show()