import csv
import math
import numpy as np


# ==========================================================
# ONLY CHANGE THIS
# ==========================================================

position = "P7"

folder = "20261013"

filename = f"{folder}/1013{position}MD.csv"


# ==========================================================
# READ CSV
# ==========================================================

with open(filename, "r", encoding="cp1252") as file:
    rows = list(csv.reader(file))


# ==========================================================
# READ PUPIL INFORMATION
# ==========================================================

pupil_diameter_x = float(rows[32][1])
pupil_diameter_y = float(rows[33][1])

pupil_center_x = float(rows[34][1])
pupil_center_y = float(rows[35][1])

Rx = pupil_diameter_x / 2
Ry = pupil_diameter_y / 2


# ==========================================================
# READ THORLABS Z1-Z15 COEFFICIENTS
#
# Excel rows 45-59
# Python rows 44-58
# ==========================================================

thorlabs = np.array([
    float(row[3])
    for row in rows[44:59]
])


# ==========================================================
# READ WAVEFRONT GRID
# ==========================================================

# X coordinates
x_values = np.array([
    float(value)
    for value in rows[100][1:10]
])


# Y coordinates + wavefront values
y_values = []
wavefront = []


for row in rows[101:110]:

    y_values.append(
        float(row[0])
    )

    current_row = []

    for value in row[1:10]:

        if value.strip() == "NaN":
            current_row.append(np.nan)

        else:
            current_row.append(
                float(value)
            )

    wavefront.append(current_row)


y_values = np.array(y_values)
wavefront = np.array(wavefront)


# ==========================================================
# CREATE X-Y GRID
# ==========================================================

X, Y = np.meshgrid(
    x_values,
    y_values
)


# ==========================================================
# KEEP VALID WAVEFRONT VALUES
# ==========================================================

valid = np.isfinite(wavefront)

x = X[valid]
y = Y[valid]
W = wavefront[valid]


# ==========================================================
# NORMALIZE COORDINATES USING ACTUAL PUPIL
# ==========================================================

xn = (x - pupil_center_x) / Rx
yn = (y - pupil_center_y) / Ry

rho = np.sqrt(
    xn**2 + yn**2
)

phi = np.arctan2(
    -yn,
    xn
)


# Only use points inside pupil
inside = rho <= 1

rho = rho[inside]
phi = phi[inside]
W = W[inside]


# ==========================================================
# RADIAL ZERNIKE POLYNOMIAL
# ==========================================================

def radial_zernike(n, m, rho):

    m = abs(m)

    R = np.zeros_like(rho)

    for k in range((n - m) // 2 + 1):

        coefficient = (
            (-1)**k
            * math.factorial(n - k)
            /
            (
                math.factorial(k)
                * math.factorial((n + m)//2 - k)
                * math.factorial((n - m)//2 - k)
            )
        )

        R += (
            coefficient
            * rho**(n - 2*k)
        )

    return R


# ==========================================================
# THORLABS / ANSI NORMALIZED ZERNIKE
# ==========================================================

def zernike(n, m, rho, phi):

    R = radial_zernike(
        n,
        m,
        rho
    )

    # m = 0
    if m == 0:

        return (
            np.sqrt(n + 1)
            * R
        )

    # m != 0
    norm = np.sqrt(
        2 * (n + 1)
    )

    if m > 0:

        return (
            norm
            * R
            * np.cos(m * phi)
        )

    else:

        return (
            norm
            * R
            * np.sin(abs(m) * phi)
        )


# ==========================================================
# THORLABS Z1-Z15
# ==========================================================

modes = [
    (0,  0),   # Z1  Piston
    (1, -1),   # Z2  Y Tilt
    (1,  1),   # Z3  X Tilt
    (2, -2),   # Z4  Oblique Astigmatism
    (2,  0),   # Z5  Defocus
    (2,  2),   # Z6
    (3, -3),   # Z7
    (3, -1),   # Z8
    (3,  1),   # Z9
    (3,  3),   # Z10
    (4, -4),   # Z11
    (4, -2),   # Z12
    (4,  0),   # Z13
    (4,  2),   # Z14
    (4,  4)    # Z15
]


# ==========================================================
# CALCULATE ZERNIKE VALUES AT EACH WAVEFRONT POINT
# ==========================================================

columns = []

for n, m in modes:

    columns.append(
        zernike(
            n,
            m,
            rho,
            phi
        )
    )


Z = np.column_stack(columns)


# ==========================================================
# LEAST-SQUARES FIT
#
# W ≈ C1*Z1 + C2*Z2 + ... + C15*Z15
# ==========================================================

ours, _, _, _ = np.linalg.lstsq(
    Z,
    W,
    rcond=None
)


# ==========================================================
# PRINT ALL 15
# ==========================================================

print()
print(f"MD-BASED ZERNIKE FIT FOR {position}")
print("------------------------------------------------")
print("Mode        Ours        Thorlabs       Difference")
print("------------------------------------------------")


for i in range(15):

    difference = ours[i] - thorlabs[i]

    print(
        f"Z{i+1:<3}       "
        f"{ours[i]:>8.4f}     "
        f"{thorlabs[i]:>8.4f}     "
        f"{difference:>8.4f}"
    )


# ==========================================================
# PRINT Z2-Z5 ONLY
# ==========================================================

print()
print("Z2-Z5 OF INTEREST")
print("------------------------------------------------")

print(
    f"Z2 - Y Tilt:              "
    f"Ours = {ours[1]:.4f}   "
    f"Thorlabs = {thorlabs[1]:.4f}"
)

print(
    f"Z3 - X Tilt:              "
    f"Ours = {ours[2]:.4f}   "
    f"Thorlabs = {thorlabs[2]:.4f}"
)

print(
    f"Z4 - Oblique Astigmatism: "
    f"Ours = {ours[3]:.4f}   "
    f"Thorlabs = {thorlabs[3]:.4f}"
)

print(
    f"Z5 - Defocus:             "
    f"Ours = {ours[4]:.4f}   "
    f"Thorlabs = {thorlabs[4]:.4f}"
)