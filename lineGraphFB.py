import csv
import matplotlib.pyplot as plt


# Files in order from original position to furthest backward position
files = [
    "20261013/1013P0MD.csv",
    "20261013/1013P1MD.csv",
    "20261013/1013P2MD.csv",
    "20261013/1013P3MD.csv",
    "20261013/1013P4MD.csv",
    "20261013/1013P5MD.csv",
    "20261013/1013P6MD.csv",
    "20261013/1013P7MD.csv",
    "20261013/1013P8MD.csv",
    "20261013/1013P9MD.csv"
]


# P0 = original calibrated position
positions = list(range(10))


def get_coefficients(filename):

    with open(filename, "r", encoding="cp1252") as file:
        rows = list(csv.reader(file))

    # Thorlabs Zernike modes
    z2 = float(rows[45][3])   # Z2 = Tip Y
    z3 = float(rows[46][3])   # Z3 = Tilt X
    z4 = float(rows[47][3])   # Z4 = Astigmatism (n=2, m=-2)
    z5 = float(rows[48][3])   # Z5 = Defocus

    return z2, z3, z4, z5


# Empty lists for each Zernike coefficient
z2_values = []
z3_values = []
z4_values = []
z5_values = []


# Read each file
for file in files:

    z2, z3, z4, z5 = get_coefficients(file)

    z2_values.append(z2)
    z3_values.append(z3)
    z4_values.append(z4)
    z5_values.append(z5)


# Print values so we can inspect them
print("Positions:", positions)
print("Z2 - Tip Y:", z2_values)
print("Z3 - Tilt X:", z3_values)
print("Z4 - Astigmatism", z4_values)
print("Z5 - Defocus:", z5_values)


# Plot
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
    label="Z4 - Astigmatism"
)

plt.plot(
    positions,
    z5_values,
    marker="o",
    label="Z5 - Defocus"
)


# Zero reference line
plt.axhline(0, linewidth=1)

# Labels
plt.xlabel("Backward Shift Position")
plt.ylabel("Zernike Coefficient (µm)")
plt.title("Zernike Coefficients vs Backward Position")

# Make x-axis say P0, P1, ..., P9
plt.xticks(
    positions,
    ["P0", "P1", "P2", "P3", "P4",
     "P5", "P6", "P7", "P8", "P9"]
)

plt.legend()
plt.tight_layout()
plt.show()