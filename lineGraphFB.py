import csv
import matplotlib.pyplot as plt

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


positions = [-2, -1, 0, 1, 2, 3, 4, 5]


def get_coefficients(filename):

    with open(filename, "r", encoding="cp1252") as file:
        rows = list(csv.reader(file))

    # Excel Index 2 = real-world Zernike 1
    zernike1 = float(rows[45][3])

    # Excel Index 3 = real-world Zernike 2
    zernike2 = float(rows[46][3])

    # Excel Index 4 = real-world Zernike 3
    zernike3 = float(rows[47][3])

    # Excel Index 5 = real-world Zernike 4
    zernike4 = float(rows[48][3])

    return zernike1, zernike2, zernike3, zernike4

z1 = []
z2 = []
z3 = []
z4 = []


for file in files:

    c1, c2, c3, c4 = get_coefficients(file)

    z1.append(c1)
    z2.append(c2)
    z3.append(c3)
    z4.append(c4)



print("Positions:", positions)
print("Zernike 1:", z1)
print("Zernike 2:", z2)
print("Zernike 3:", z3)
print("Zernike 4:", z4)



plt.figure(figsize=(10, 6))

plt.plot(
    positions,
    z1,
    marker="o",
    label="Zernike 1"
)

plt.plot(
    positions,
    z2,
    marker="o",
    label="Zernike 2"
)

plt.plot(
    positions,
    z3,
    marker="o",
    label="Zernike 3"
)

plt.plot(
    positions,
    z4,
    marker="o",
    label="Zernike 4"
)



plt.axhline(0, linewidth=1)


plt.xlabel("Position")
plt.ylabel("Zernike Coefficient (µm)")
plt.title("Zernike Coefficients vs Forward/Backward Position")


plt.xticks(positions)



plt.legend()

plt.tight_layout()
plt.show()