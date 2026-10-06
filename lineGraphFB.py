import csv
import matplotlib.pyplot as plt


# Files in position order from -2 to +5
files = [
    "P-2_MD.csv",
    "P-1_MD.csv",
    "P0_MD.csv",
    "P1_MD.csv",
    "P2_MD.csv",
    "P3_MD.csv",
    "P4_MD.csv",
    "P5_MD.csv"
]


# Corresponding positions
positions = [-2, -1, 0, 1, 2, 3, 4, 5]


# Get Zernike coefficients from each file
def get_coefficients(filename):

    with open(filename, "r", encoding="cp1252") as file:
        rows = list(csv.reader(file))

    # Thorlabs Z2 = Tip Y
    z2 = float(rows[45][3])

    # Thorlabs Z3 = Tilt X
    z3 = float(rows[46][3])

    # Thorlabs Z4 = Astigmatism 45°
    z4 = float(rows[47][3])

    # Thorlabs Z5 = Defocus
    z5 = float(rows[48][3])

    return z2, z3, z4, z5


# Empty lists
z2_values = []
z3_values = []
z4_values = []
z5_values = []


# Extract coefficients
for file in files:

    z2, z3, z4, z5 = get_coefficients(file)

    z2_values.append(z2)
    z3_values.append(z3)
    z4_values.append(z4)
    z5_values.append(z5)


# Print values
print("Positions:", positions)
print("Z2 - Tip Y:", z2_values)
print("Z3 - Tilt X:", z3_values)
print("Z4 - Astigmatism 45°:", z4_values)
print("Z5 - Defocus:", z5_values)


# Make line plot
plt.figure(figsize=(10, 6))

plt.plot(
    positions,
    z2_values,
    marker="o",
    label="Z2 - Tip Y"
)

plt.plot(
    positions,
    z3_values,
    marker="o",
    label="Z3 - Tilt X"
)

plt.plot(
    positions,
    z4_values,
    marker="o",
    label="Z4 - Astigmatism 45°"
)

plt.plot(
    positions,
    z5_values,
    marker="o",
    label="Z5 - Defocus"
)


plt.axhline(0, linewidth=1)

plt.xlabel("Position")
plt.ylabel("Zernike Coefficient (µm)")
plt.title("Zernike Coefficients vs Forward/Backward Position")

plt.xticks(positions)

plt.legend()

plt.tight_layout()
plt.show()