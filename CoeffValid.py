import csv
import math
import numpy as np


# ==========================================================
# ONLY CHANGE THIS
# ==========================================================

position = "P8"

folder = "20261013"

reference_file = f"{folder}/20261013ref(Sheet1).csv"
centroid_file = f"{folder}/1013{position}CD.csv"
measurement_file = f"{folder}/1013{position}MD.csv"


# ==========================================================
# HARDWARE CONSTANTS
# WFS300-14AR / MLA300M-14AR
# ==========================================================

pixel_size_mm = 0.00465      # DCU224M: 4.65 µm pixel
focal_length_mm = 14.2       # MLA300M-14AR effective focal length


# ==========================================================
# READ CSV
# ==========================================================

def read_csv(filename):

    with open(filename, "r", encoding="cp1252") as file:
        return list(csv.reader(file))


# ==========================================================
# READ EXACT USER CALIBRATION REFERENCE
# ==========================================================

ref_rows = read_csv(reference_file)

ref_x = np.full((9, 9), np.nan)
ref_y = np.full((9, 9), np.nan)


for r in range(9):

    for c in range(9):

        x = float(ref_rows[r][2*c])
        y = float(ref_rows[r][2*c + 1])

        # 0,0 means there is no reference spot
        if x == 0 and y == 0:
            continue

        ref_x[r, c] = x
        ref_y[r, c] = y


# ==========================================================
# READ CURRENT CENTROIDS
# ==========================================================

cd_rows = read_csv(centroid_file)

meas_x = np.full((9, 9), np.nan)
meas_y = np.full((9, 9), np.nan)


# Excel rows 42-50
# Python rows 41-49

for r in range(9):

    row = cd_rows[41 + r]

    for c in range(9):

        try:

            meas_x[r, c] = float(row[2*c])
            meas_y[r, c] = float(row[2*c + 1])

        except:
            pass


# ==========================================================
# READ MD FILE
#
# IMPORTANT:
# We are NOT using its wavefront values to calculate
# coefficients.
#
# We only use:
#   1. Pupil settings
#   2. Physical x/y coordinates
#   3. Valid-spot mask
#   4. Thorlabs coefficients for comparison
# ==========================================================

md_rows = read_csv(measurement_file)


# ==========================================================
# PUPIL USED BY THORLABS FOR THIS MD MEASUREMENT
# ==========================================================

pupil_diameter_x = float(md_rows[32][1])
pupil_diameter_y = float(md_rows[33][1])

pupil_center_x = float(md_rows[34][1])
pupil_center_y = float(md_rows[35][1])


Rx = pupil_diameter_x / 2
Ry = pupil_diameter_y / 2


# ==========================================================
# THORLABS COEFFICIENTS
# ==========================================================

thorlabs_z2 = float(md_rows[45][3])
thorlabs_z3 = float(md_rows[46][3])
thorlabs_z4 = float(md_rows[47][3])
thorlabs_z5 = float(md_rows[48][3])


# ==========================================================
# READ PHYSICAL SPOT COORDINATES DIRECTLY FROM MD FILE
#
# Do NOT manually invent ±1.2 grid orientation anymore.
# ==========================================================

x_coordinates = np.array([
    float(value)
    for value in md_rows[100][1:10]
])


y_coordinates = np.array([
    float(md_rows[101 + r][0])
    for r in range(9)
])


X, Y = np.meshgrid(
    x_coordinates,
    y_coordinates
)


# ==========================================================
# READ MD VALID-SPOT MASK
#
# We do NOT use the wavefront numbers themselves.
# NaN simply tells us that the spot was not valid in the
# acquisition that produced the Thorlabs Zernike values.
# ==========================================================

md_valid = np.zeros((9, 9), dtype=bool)


for r in range(9):

    for c in range(9):

        try:

            value = float(
                md_rows[101 + r][1 + c]
            )

            md_valid[r, c] = np.isfinite(value)

        except:

            md_valid[r, c] = False


# ==========================================================
# SPOT DEVIATIONS
#
# The saved centroid data use camera pixel coordinates.
#
# WFS physical coordinate convention is viewed toward
# the source, while the raw image coordinates use the
# camera array convention.
#
# Therefore use reference - measurement here.
# ==========================================================

delta_x_pixels = ref_x - meas_x
delta_y_pixels = ref_y - meas_y


# ==========================================================
# CONVERT SPOT DEVIATIONS TO WAVEFRONT SLOPES
#
# dW/dx = delta_x / f
# dW/dy = delta_y / f
#
# ×1000 converts from mm/mm to µm/mm
# ==========================================================

slope_x = (
    1000
    * delta_x_pixels
    * pixel_size_mm
    / focal_length_mm
)


slope_y = (
    1000
    * delta_y_pixels
    * pixel_size_mm
    / focal_length_mm
)


# ==========================================================
# VALID CURRENT / REFERENCE CENTROIDS
# ==========================================================

centroid_valid = (
    np.isfinite(ref_x)
    &
    np.isfinite(ref_y)
    &
    np.isfinite(meas_x)
    &
    np.isfinite(meas_y)
)


# ==========================================================
# THORLABS PUPIL MASK
# ==========================================================

inside_pupil = (
    ((X - pupil_center_x) / Rx)**2
    +
    ((Y - pupil_center_y) / Ry)**2
    <= 1
)


# Use only spots that:
#
# 1. exist in exact reference
# 2. exist in our CD measurement
# 3. were valid during MD measurement
# 4. lie inside Thorlabs pupil
# ==========================================================

use = (
    centroid_valid
    &
    md_valid
    &
    inside_pupil
)


x = X[use]
y = Y[use]

sx = slope_x[use]
sy = slope_y[use]


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

        R += coefficient * rho**(n - 2*k)

    return R


# ==========================================================
# THORLABS NORMALIZED ZERNIKE
# ==========================================================

def zernike(n, m, x, y):

    xn = (x - pupil_center_x) / Rx
    yn = (y - pupil_center_y) / Ry


    rho = np.sqrt(
        xn**2 + yn**2
    )


    phi = np.arctan2(
        yn,
        xn
    )


    R = radial_zernike(
        n,
        m,
        rho
    )


    if m == 0:

        return (
            np.sqrt(n + 1)
            * R
        )


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
# FIT Z2-Z15
#
# Z1 piston has zero derivative, so slope data cannot
# determine it.
# ==========================================================

modes = [
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
# CALCULATE ZERNIKE DERIVATIVES
# ==========================================================

h = 0.00001

columns = []


for n, m in modes:

    dZdx = (
        zernike(n, m, x + h, y)
        -
        zernike(n, m, x - h, y)
    ) / (2*h)


    dZdy = (
        zernike(n, m, x, y + h)
        -
        zernike(n, m, x, y - h)
    ) / (2*h)


    columns.append(
        np.concatenate([
            dZdx,
            dZdy
        ])
    )


# ==========================================================
# LEAST-SQUARES FIT
# ==========================================================

A = np.column_stack(columns)


measured_slopes = np.concatenate([
    sx,
    sy
])


coefficients, _, _, _ = np.linalg.lstsq(
    A,
    measured_slopes,
    rcond=None
)


# ==========================================================
# Z2-Z5
# ==========================================================

our_z2 = coefficients[0]
our_z3 = coefficients[1]
our_z4 = coefficients[2]
our_z5 = coefficients[3]


# ==========================================================
# OUTPUT
# ==========================================================

print()
print(f"CENTROID-BASED ZERNIKE FIT FOR {position}")
print("----------------------------------------------------")

print(f"CD acquisition time: {cd_rows[6][1].strip()}")
print(f"MD acquisition time: {md_rows[6][1].strip()}")

print()

print(f"Spots used: {len(x)}")

print()

print(
    f"Z2 - Y Tilt:              "
    f"Ours = {our_z2:.4f}   "
    f"Thorlabs = {thorlabs_z2:.4f}"
)

print(
    f"Z3 - X Tilt:              "
    f"Ours = {our_z3:.4f}   "
    f"Thorlabs = {thorlabs_z3:.4f}"
)

print(
    f"Z4 - Oblique Astigmatism: "
    f"Ours = {our_z4:.4f}   "
    f"Thorlabs = {thorlabs_z4:.4f}"
)

print(
    f"Z5 - Defocus:             "
    f"Ours = {our_z5:.4f}   "
    f"Thorlabs = {thorlabs_z5:.4f}"
)