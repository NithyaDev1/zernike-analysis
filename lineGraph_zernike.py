import csv
import matplotlib.pyplot as plt


files = [
    "MD_pos0.csv",
    "MD_pos1.csv",
    "MD_pos2.csv",
    "MD_pos3.csv",
    "MD_pos4.csv",
    "MD_pos5.csv",
    "MD_pos6.csv",
    "MD_pos7.csv",
    "MD_pos8.csv",
    "MD_pos9.csv",
    "MD_pos10.csv"
]


def get_coefficients(filename):

    with open(filename, "r", encoding="cp1252") as file:
        rows = list(csv.reader(file))

    coefficient2 = float(rows[45][3])
    coefficient3 = float(rows[46][3])
    coefficient4 = float(rows[47][3])
    coefficient5 = float(rows[48][3])

    return coefficient2, coefficient3, coefficient4, coefficient5


positions = list(range(11))

coeff2 = []
coeff3 = []
coeff4 = []
coeff5 = []


for file in files:

    c2, c3, c4, c5 = get_coefficients(file)

    coeff2.append(c2)
    coeff3.append(c3)
    coeff4.append(c4)
    coeff5.append(c5)


plt.figure(figsize=(10, 6))

plt.plot(
    positions,
    coeff2,
    marker="o",
    label="Zernike Index 2",
    color="blue"
)

plt.plot(
    positions,
    coeff3,
    marker="o",
    label="Zernike Index 3",
    color="red"
)

plt.plot(
    positions,
    coeff4,
    marker="o",
    label="Zernike Index 4",
    color="green"
)

plt.plot(
    positions,
    coeff5,
    marker="o",
    label="Zernike Index 5",
    color="purple"
)

plt.axhline(0, linewidth=1)

plt.xlabel("Position")
plt.ylabel("Zernike Coefficient (µm)")
plt.title("Zernike Coefficients vs Position")

plt.xticks(positions)

plt.legend()

plt.tight_layout()
plt.show()