##################################################################
# cylinder.m: Channel flow past a cylinderical        
#             obstacle, using a LB method            
##################################################################
# Lattice Boltzmann sample in Matlab
# Copyright (C) 2006-2008 Jonas Latt
# Address: EPFL, 1015 Lausanne, Switzerland
# E-mail: jonas@lbmethod.org
# Get the most recent version of this file on LBMethod.org:
#   http://www.lbmethod.org/_media/numerics:cylinder.m
#
# Original implementaion of Zou/He boundary condition by
# Adriano Sciacovelli (see example "cavity.m")
##################################################################
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public 
# License along with this program; if not, write to the Free 
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
##################################################################

import numpy as np
import matplotlib.pyplot as plt
import scipy.io

# GENERAL FLOW CONSTANTS
lx = 400                        # number of cells in x-direction
ly = 100                        # number of cells in y-direction

lx = 2 * 50                     # number of cells in x-direction
ly = 2 * 25                     # number of cells in y-direction

obst_x = lx / 5 + 1             # position of the cylinder (exact
obst_y = ly / 2 + 3             # y-symmetry is avoided)
obst_r = ly / 10 + 1            # radius of the cylinder
uMax = 0.1                      # maximum velocity of Poiseuille inflow
Re = 100                        # Reynolds number
nu = uMax * 2 * obst_r / Re     # kinematic viscosity
omega = 1 / (3 * nu + 1 / 2)    # relaxation parameter
# maxT   = 2*400000             # total number of iterations
tPlot = 1 * 50                  # cycles
maxT = 5000                     # total number of iterations

# up = zeros([1, maxT])
up = np.zeros(maxT)

# D2Q9 LATTICE CONSTANTS
t = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
cx = [0, 1, 0, -1, 0, 1, -1, -1, 1]
cy = [0, 0, 1, 0, -1, 1, 1, -1, -1]
opp = [1, 4, 5, 2, 3, 8, 9, 6, 7]
# col = [2:(ly - 1)]
col = np.arange(2, ly)
inlet = 1  # position of inlet
outlet = lx  # position of outlet

# [y,x] = meshgrid(1:ly,1:lx) # get coordinate of matrix indices
y, x = np.meshgrid(np.arange(1, ly + 1), np.arange(1, lx + 1))
# x is downstream.

# obst = ...                   # Location of cylinder
#     (x-obst_x).^2 + (y-obst_y).^2 <= obst_r.^2
# obst(:,[1,ly]) = 1    # Location of top/bottom boundary
# bbRegion = find(obst) # Boolean mask for bounce-back cells

obst = (x - obst_x) ** 2 + (y - obst_y) ** 2 <= obst_r ** 2  # Mask of the location of the cylinder
obst[:, [0, ly - 1]] = True  # Set location of top/bottom boundary
# bbRegion = np.where(obst)  # Boolean mask for bounce-back cells
bbRegion = np.where(obst.T.flatten())  # Same behaviour as matlab
# np.nonzero(obst.T.flatten())[0]

# INITIAL CONDITION: Poiseuille profile at equilibrium
L = ly - 2
y_phys = y - 1.5
ux = 4 * uMax / (L * L) * (y_phys * L - y_phys * y_phys)
# uy = zeros(lx,ly)
uy = np.zeros((lx, ly))
rho = 1
fIn = np.zeros((9, lx, ly))

for i in range(9):
	cu = 3 * (cx[i] * ux + cy[i] * uy)
	fIn[i] = rho * t[i] * (1 + cu + 1 / 2 * (cu * cu) - 3 / 2 * (ux ** 2 + uy ** 2))

# MAIN LOOP (TIME CYCLES)
for cycle in range(maxT):

	# MACROSCOPIC VARIABLES
	rho = sum(fIn)
	# np.sum(fIn, axis=0) would also work
	fIn_reshape = np.reshape(fIn, (9, lx * ly), order="F")
	# Matlab is fortran style for reshape.
	# This is also a solution.
	# fIn_reshape = fIn.T.reshape((lx*ly, 9)).T
	ux = np.reshape(np.dot(cx, fIn_reshape), (1, lx, ly), order="F") / rho

	# TODO: uy is somewhat different than MATLAB. Assuming it is do to precision
	uy = np.reshape(np.dot(cy, fIn_reshape), (1, lx, ly), order="F") / rho
	# For some reason the first column as been put last compare to the matlab
	uy_last = uy[..., -1].copy()
	uy[..., 1:] = uy[..., :-1]
	uy[..., 0] = uy_last

	# MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS
	# Inlet: Poiseuille profile
	y_phys = col - 1.5

	ux[:, inlet, col] = 4 * uMax / (L * L) * (y_phys * L - y_phys * y_phys)
	uy[:, inlet, col] = 0
	rho[:, inlet, col] = 1 / (1 - ux[:, inlet, col]) * (sum(fIn[[0, 2, 4], inlet, col]) + 2 * sum(fIn[[5, 6, 7], inlet, col]))

	# Outlet: Constant pressure
	rho[:, outlet, col] = 1
	ux[:, outlet, col] = -1 + 1 / (rho[:, outlet, col]) * (sum(fIn[[0, 2, 4], outlet, col]) + 2 * sum(fIn[[1, 5, 8], outlet, col]))
	uy[:, outlet, col] = 0

	# MICROSCOPIC BOUNDARY CONDITIONS: INLET (Zou/He BC)
	fIn[1, inlet, col] = fIn[3, inlet, col] + 2 / 3 * rho[:, inlet, col] * ux[:, inlet, col]

	fIn[5, inlet, col] = (
			fIn[7, inlet, col] +
			1 / 2 * (fIn[4, inlet, col] - fIn[2, inlet, col]) +
			1 / 2 * rho[:, inlet, col] * uy[:, inlet, col] +
			1 / 6 * rho[:, inlet, col] * ux[:, inlet, col]
	)

	fIn[8, inlet, col] = (
			fIn[6, inlet, col] +
			1 / 2 * (fIn[2, inlet, col] - fIn[4, inlet, col]) -
			1 / 2 * rho[:, inlet, col] * uy[:, inlet, col] +
			1 / 6 * rho[:, inlet, col] * ux[:, inlet, col]
	)

	# MICROSCOPIC BOUNDARY CONDITIONS: OUTLET (Zou/He BC)
	fIn[3, outlet, col] = fIn[1, outlet, col] - 2 / 3 * rho[:, outlet, col] * ux[:, outlet, col]
	fIn[7, outlet, col] = (
			fIn[6, outlet, col] +
			1 / 2 * (fIn[2, outlet, col] - fIn[5, outlet, col]) -
			1 / 2 * rho[:, outlet, col] * uy[:, outlet, col] -
			1 / 6 * rho[:, outlet, col] * ux[:, outlet, col]
	)
	fIn[6, outlet, col] = (
			fIn[8, outlet, col] +
			1 / 2 * (fIn[4, outlet, col] - fIn[2, outlet, col]) +
			1 / 2 * rho[:, outlet, col] * uy[:, outlet, col] -
			1 / 6 * rho[:, outlet, col] * ux[:, outlet, col]
	)

	fEq = np.zeros((9, lx, ly))  # Initialize equilibrium distribution
	fOut = np.zeros((9, lx, ly))  # Initialize output distribution

	# COLLISION STEP
	for i in range(9):
		cu = 3 * (cx[i] * ux + cy[i] * uy)
		fEq[i, :, :] = rho * t[i] * (1 + cu + 1 / 2 * (cu * cu) - 3 / 2 * (ux ** 2 + uy ** 2))
		fOut[i, :, :] = fIn[i, :, :] - omega * (fIn[i, :, :] - fEq[i, :, :])

	# OBSTACLE (BOUNCE-BACK)
	for i in range(9):
		fOut[i, bbRegion] = fIn[opp[i], bbRegion]

	# STREAMING STEP
	for i in range(9):
		fIn[i, :, :] = np.roll(fOut[i, :, :], shift=(cx[i], cy[i]), axis=0)

	# VISUALIZATION
	if cycle % tPlot == 0:
		print(cycle)

	u = np.reshape(np.sqrt(ux ** 2 + uy ** 2), (lx, ly), order="F")
	u[bbRegion] = np.NAN
	plt.imshow(u.conj().T)
	plt.axis("equal")
	plt.axis("off")
	plt.show()

	# printf('%d of %d \n', cycle, maxT)
	vx = np.reshape(ux, (lx, ly), order="F")
	vy = np.reshape(uy, (lx, ly), order="F")
	fileName = f"Testvelocity_{cycle}.mat"
	data = {"vx": vx, "vy": vy}
	scipy.io.savemat(fileName, data)

	vx = np.reshape(ux, (lx, ly), order="F")
	vy = np.reshape(uy, (lx, ly), order="F")
	# no idea what the values in vx(..,..) should be. Just picking as smaller than the overall window size
	u = np.sqrt(vx[2 * 50, 2 * 25] ** 2 + vy[2 * 50, 2 * 25] ** 2)
	up[cycle + 2] = u

data = {
	"bbRegion":  bbRegion,
	"x":  x,
	"y":  y,
	"obst_x":  obst_x,
	"obst_y":  obst_y,
	"obst_r":  obst_r,
	"nu":  nu,
	"Re":  Re,
	"uMax":  uMax,
	"lx":  lx,
	"ly":  ly,
}
np.save(
	'TestSymVars.mat',
	data
)


def main():
	pass


if __name__ == '__main__':
	main()
