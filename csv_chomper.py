import numpy as np
import csv
from math import sqrt

def read_csv_fields(filename):
    """
    Reads the CSV file and returns:
    - u[t, x]: velocity field
    - rho[t, x]: density field from |psi|²
    - times[t]: time values
    """
    u_list = []
    rho_list = []
    times = []

    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        row_iter = iter(reader)

        for psi_row in row_iter:
            try:
                u_row = next(row_iter)
            except StopIteration:
                break

            # Parse psi_row
            psi_row = [float(v) for v in psi_row]
            time = psi_row[0]
            real = np.array(psi_row[1::2])
            imag = np.array(psi_row[2::2])
            rho = real**2 + imag**2

            # Parse u_row
            u_row = [float(v) for v in u_row]
            u = np.array(u_row[1:])  # skip time column

            times.append(time)
            rho_list.append(rho)
            u_list.append(u)

    return np.array(u_list), np.array(rho_list), np.array(times)


def c(rho, gamma):
    return sqrt(gamma-1)*rho**((gamma-1)/2)

def lambda_plus(rho, u, gamma):
    return u + 2/(gamma-1) * c(rho, gamma)

def lambda_minus(rho, u, gamma):
    return u - 2/(gamma-1) * c(rho, gamma)

#supposes that the bulk density is one (so speed of sound is 1)
def lambda_zero(gamma):
    return 2.0/(gamma -1.0) * c(1., gamma)

import numpy as np
def lambda_moins_theo(L, gamma, t, x):
    l_0 = lambda_zero(gamma)
    x_moins = L - c(1., gamma) * t - 1e-5
    x_plus  = L + l_0 * t + 1e-5  # Avoid edge-case singularities

    # Linear interpolation formula:
    slope = 2 * l_0 / (x_plus - x_moins)
    midpoint = (x_plus + x_moins) / 2
    linear_interp = slope * (x - midpoint)

    result = np.where(
        x < x_moins,
        -l_0,
        np.where(
            x > x_plus,
            l_0,
            linear_interp
        )
    )
    return result

def lambda_plus_theo(L, gamma, t, x):
    l_0 = lambda_zero(gamma)
    x_plus = -L + c(1., gamma) * t + 1e-5
    x_moins = -L - lambda_zero(gamma) * t - 1e-5
    t_coll = L/c(1., gamma)
  # Linear interpolation formula:
    slope = 2 * l_0 / (x_plus - x_moins)
    midpoint = (x_plus + x_moins) / 2
    linear_interp = slope * (x - midpoint)

    result = np.where(
        x < x_moins,
        -l_0,
        np.where(
            x > x_plus,
            l_0,
            linear_interp
        )
    )
    return result

def lambda_to_urho(gamma, lambda_plus, lambda_moins):
    u = 1/2 * (lambda_plus + lambda_moins)
    c = (gamma - 1)/4 * (lambda_plus - lambda_moins)
    rho = (1/sqrt(gamma-1.) * c)**(2/(gamma-1))
    return rho, u


