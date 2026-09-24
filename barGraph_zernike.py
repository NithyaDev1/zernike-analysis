import csv
import matplotlib.pyplot as plt
import numpy as np


files = [
    "test_measurement.csv",
    "position2_MD.csv",
    "position3_MD.csv",
    "position4_MD.csv",
    "position5_MD.csv"
]

labels = [
    "test_measurement",
    "Position 2",
    "Position 3",
    "Position 4",
    "Position 5"
]


def get_coefficients(filename):

    with open(filename, "r", encoding="cp1252") as file:
        rows = list(csv.reader(file))

    # Excel rows 45 through 59
    zernike_rows = rows[44:59]

    # First column = Index
    indices = [int(row[0]) for row in zernike_rows]

    # Fourth column = Coefficient
    coefficients = [float(row[3]) for row in zernike_rows]

    return indices, coefficients


# Spacing between Zernike groups
x = np.arange(15) * 3

# Width of each individual bar
width = 0.5

plt.figure(figsize=(16, 7))


for i in range(len(files)):

    indices, coefficients = get_coefficients(files[i])

    position = x + (i - 2) * width

    plt.bar(
        position,
        coefficients,
        width=width,
        label=labels[i]
    )


plt.axhline(0, linewidth=1)

plt.xlabel("Zernike Index")
plt.ylabel("Coefficient (µm)")
plt.title("Zernike Coefficient Comparison")

plt.xticks(x, indices)

plt.legend()
plt.tight_layout()

plt.show()