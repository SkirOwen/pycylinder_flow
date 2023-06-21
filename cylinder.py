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

# Copyright (C) 2023  Owen Allemang
#
# This file is part of PyCylinder_flow.
# PyCylinder_flow is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# PyCylinder_flow is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Foobar.
# If not, see <https://www.gnu.org/licenses/>.

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from PIL import Image
from tqdm import tqdm

import cmocean as cm

from utils import clear_dir, plot_dir, velocity_dir


def simulate(lx: int = 400, ly: int = 100, max_t: int = 800_000):
	obst_x = lx / 5 + 1             # position of the cylinder (exact
	obst_y = ly / 2 + 3             # y-symmetry is avoided)
	obst_r = ly / 10 + 1            # radius of the cylinder

	uMax = 0.1                      # maximum velocity of Poiseuille inflow
	Re = 100                        # Reynolds number
	nu = uMax * 2 * obst_r / Re     # kinematic viscosity
	omega = 1 / (3 * nu + 1 / 2)    # relaxation parameter

	tPlot = 1 * 50                  # cycles

	# up = np.zeros(max_t)

	# D2Q9 LATTICE CONSTANTS
	t = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])
	cx = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
	cy = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])
	opp = [0, 3, 4, 1, 2, 7, 8, 5, 6]
	# col = [2:(ly - 1)]
	col = np.arange(1, ly-1)
	inlet = 0               # position of inlet
	outlet = lx - 1         # position of outlet

	# x is downstream
	y, x = np.meshgrid(np.arange(1, ly + 1), np.arange(1, lx + 1))

	# Mask of the location of the cylinder
	obst = (x - obst_x) ** 2 + (y - obst_y) ** 2 <= obst_r ** 2
	# Set location of top/bottom boundary
	obst[:, [0, ly - 1]] = True
	# Boolean mask for bounce-back cells
	bbRegion = np.where(obst)

	# INITIAL CONDITION: Poiseuille profile at equilibrium
	L = ly - 2
	y_phys = y - 1.5
	ux = 4 * uMax / (L * L) * (y_phys * L - y_phys * y_phys)
	uy = np.zeros((lx, ly))
	rho = 1
	fIn = np.zeros((9, lx, ly))

	for i in range(9):
		cu = 3 * (cx[i] * ux + cy[i] * uy)
		fIn[i] = rho * t[i] * (1 + cu + 1/2 * (cu * cu) - 3/2 * (ux ** 2 + uy ** 2))

	# MAIN LOOP (TIME CYCLES)
	for cycle in tqdm(range(max_t), desc="Lattice Boltzmann Simulation", unit="step"):

		# MACROSCOPIC VARIABLES
		rho = np.sum(fIn, axis=0)
		fIn_reshape = np.reshape(fIn, (9, lx * ly), order="F")
		# Matlab is fortran style for reshape.
		ux = np.reshape(np.dot(cx, fIn_reshape), (lx, ly), order="F") / rho
		# uy is somewhat different from MATLAB. Assuming it is due to precision
		uy = np.reshape(np.dot(cy, fIn_reshape), (lx, ly), order="F") / rho
		# For some reason the first column as been put last compare to the matlab

		# MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS
		# Inlet: Poiseuille profile
		y_phys = col - 0.5

		ux[inlet, col] = 4 * uMax / (L * L) * (y_phys * L - y_phys * y_phys)
		uy[inlet, col] = 0
		rho[inlet, col] = (
			1 / (1 - ux[inlet, col]) *
			(
				np.sum(fIn[np.ix_([0, 2, 4], [inlet], col)], axis=0) +
				2 * np.sum(fIn[np.ix_([3, 6, 7], [inlet], col)], axis=0)
			)
		)

		# Outlet: Constant pressure
		rho[outlet, col] = 1
		ux[outlet, col] = (
			-1 + 1 / (rho[outlet, col]) *
			(
				np.sum(fIn[np.ix_([0, 2, 4], [outlet], col)], axis=0) +
				2 * np.sum(fIn[np.ix_([1, 5, 8], [outlet], col)], axis=0)
			)
		)
		uy[outlet, col] = 0

		# MICROSCOPIC BOUNDARY CONDITIONS:
		# INLET (Zou/He BC)
		fIn[1, inlet, col] = fIn[3, inlet, col] + 2 / 3 * rho[inlet, col] * ux[inlet, col]
		fIn[5, inlet, col] = (
				fIn[7, inlet, col] +
				1 / 2 * (fIn[4, inlet, col] - fIn[2, inlet, col]) +
				1 / 2 * rho[inlet, col] * uy[inlet, col] +
				1 / 6 * rho[inlet, col] * ux[inlet, col]
		)
		fIn[8, inlet, col] = (
				fIn[6, inlet, col] +
				1 / 2 * (fIn[2, inlet, col] - fIn[4, inlet, col]) -
				1 / 2 * rho[inlet, col] * uy[inlet, col] +
				1 / 6 * rho[inlet, col] * ux[inlet, col]
		)

		# OUTLET (Zou/He BC)
		fIn[3, outlet, col] = fIn[1, outlet, col] - 2 / 3 * rho[outlet, col] * ux[outlet, col]
		fIn[7, outlet, col] = (
				fIn[5, outlet, col] +
				1 / 2 * (fIn[2, outlet, col] - fIn[4, outlet, col]) -
				1 / 2 * rho[outlet, col] * uy[outlet, col] -
				1 / 6 * rho[outlet, col] * ux[outlet, col]
		)
		fIn[6, outlet, col] = (
				fIn[8, outlet, col] +
				1 / 2 * (fIn[4, outlet, col] - fIn[2, outlet, col]) +
				1 / 2 * rho[outlet, col] * uy[outlet, col] -
				1 / 6 * rho[outlet, col] * ux[outlet, col]
		)

		# Initialize equilibrium distribution
		fEq = np.zeros((9, lx, ly))
		# Initialize output distribution
		fOut = np.zeros((9, lx, ly))

		# COLLISION STEP
		for i in range(9):
			cu = 3 * (cx[i] * ux + cy[i] * uy)
			fEq[i, :, :] = rho * t[i] * (1 + cu + 1/2 * (cu * cu) - 3/2 * (ux ** 2 + uy ** 2))
			fOut[i, :, :] = fIn[i, :, :] - omega * (fIn[i, :, :] - fEq[i, :, :])

		# OBSTACLE (BOUNCE-BACK)
		for i in range(9):
			fOut[i, bbRegion[0], bbRegion[1]] = fIn[opp[i], bbRegion[0], bbRegion[1]]

		# STREAMING STEP
		for i in range(9):
			fIn[i, :, :] = np.roll(fOut[i, :, :], shift=(cx[i], cy[i]), axis=(0, 1))

		save_flow_png(ux, uy, lx, ly, bbRegion, cycle, max_t)

		save_velocity(ux, uy, lx, ly, cycle)


def plot(
		ux: np.ndarray,
		uy: np.ndarray,
		lx: int,
		ly: int,
		bbRegion: np.ndarray,
		cycle: int,
		maxT: int,
		cmap=cm.cm.speed,
		save: bool = False
) -> None:
	u = np.reshape(np.sqrt(ux ** 2 + uy ** 2), (lx, ly), order="F")
	u[bbRegion[0], bbRegion[1]] = np.NAN
	plt.imshow(u.conj().T, cmap=cmap)
	plt.axis("equal")
	plt.axis("off")
	# print("A")
	plt.show()
	if save:
		plt.savefig(
		os.path.join(
			plot_dir(),
			f"{str(cycle).zfill(len(str(maxT)))}.png"
		))


def save_flow_png(
		ux: np.ndarray,
		uy: np.ndarray,
		lx: int,
		ly: int,
		bbRegion: np.ndarray,
		cycle: int,
		maxT: int,
		cmap=cm.cm.ice
) -> None:
	u = np.reshape(np.sqrt(ux ** 2 + uy ** 2), (lx, ly), order="F")
	u[bbRegion[0], bbRegion[1]] = np.NAN

	# Scale the velocity field to 0-255 range
	u_scaled = (u - np.nanmin(u)) / (np.nanmax(u) - np.nanmin(u))
	u_scaled = cmap(u_scaled.T) * 255
	u_scaled = np.uint8(u_scaled)

	# Create PIL Image object
	img = Image.fromarray(u_scaled)

	scale = 10
	new_width = int(img.width * scale)
	new_height = int(img.height * scale)

	# Resize the image
	resized_image = img.resize((new_width, new_height), resample=Image.NEAREST)

	# Save the image
	resized_image.save(
		os.path.join(
			plot_dir(),
			f"{str(cycle).zfill(len(str(maxT)))}.png"
		)
	)


def save_velocity(ux: np.ndarray, uy: np.ndarray, lx: int, ly: int, cycle: int) -> None:
	# printf('%d of %d \n', cycle, maxT)
	vx = np.reshape(ux, (lx, ly), order="F")
	vy = np.reshape(uy, (lx, ly), order="F")
	filepath = os.path.join(
		velocity_dir(),
		f"Testvelocity_{cycle}.mat"
	)

	data = {"vx": vx, "vy": vy}
	scipy.io.savemat(filepath, data)


def save_sym_mat(
		ux: np.ndarray,
		uy: np.ndarray,
		lx: int,
		ly: int,
		bbRegion: np.ndarray,
		x: np.array,
		y: np.array,
		obst_x: float,
		obst_y: float,
		obst_r: float,
		nu: float,
		Re: float,
		uMax: float,
) -> None:
	vx = np.reshape(ux, (lx, ly), order="F")
	vy = np.reshape(uy, (lx, ly), order="F")
	# # no idea what the values in vx(..,..) should be. Just picking as smaller than the overall window size
	# TODO: change this mess of indices
	u = np.sqrt(vx[(2 * 50) - 1, (2 * 25) - 1] ** 2 + vy[(2 * 50) - 1, (2 * 25) - 1] ** 2)
	# up[cycle] = u
#
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
		# ans
		# j
		# DATA
	}
	scipy.io.savemat(
		'TestSymVars.mat',
		data
	)


def main():
	clear_dir(target="./plots/", extension="png")
	simulate(max_t=400, lx=100, ly=50)


if __name__ == '__main__':
	main()
