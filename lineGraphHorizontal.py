import csv
import matplotlib.pyplot as plt


files = [
    "left5_MD.csv",
    "left4_MD.csv",
    "left3_MD.csv",
    "left2_MD.csv",
    "left1_MD.csv",
    "OG_MD.csv",
    "right1_MD.csv",
    "right2_MD.csv",
    "right3_MD.csv",
    "right4_MD.csv",
    "right5_MD.csv"
]



positions = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]



def get_coefficients(filename):

    with open(filename, "r", encoding="cp1252") as file:
        rows = list(csv.reader(file))

    coefficient2 = float(rows[45][3])
    coefficient3 = float(rows[46][3])
    coefficient4 = float(rows[47][3])
    coefficient5 = float(rows[48][3])

    return coefficient2, coefficient3, coefficient4, coefficient5



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



print("Positions:", positions)
print("Coefficient 2:", coeff2)
print("Coefficient 3:", coeff3)
print("Coefficient 4:", coeff4)
print("Coefficient 5:", coeff5)



plt.figure(figsize=(10, 6))

plt.plot(
    positions,
    coeff2,
    marker="o",
    label="Zernike Index 1"
)

plt.plot(
    positions,
    coeff3,
    marker="o",
    label="Zernike Index 2"
)

plt.plot(
    positions,
    coeff4,
    marker="o",
    label="Zernike Index 3"
)

plt.plot(
    positions,
    coeff5,
    marker="o",
    label="Zernike Index 4"
)


plt.axhline(0, linewidth=1)



plt.xlabel("Position")
plt.ylabel("Zernike Coefficient (µm)")
plt.title("Zernike Coefficients vs Position")


plt.xticks(positions)


plt.legend()

plt.tight_layout()

plt.show()