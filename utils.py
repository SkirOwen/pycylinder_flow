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

import glob
import os


def guarantee_existence(path: str) -> str:
	"""Function to guarantee the existence of a path, and returns its absolute path.

	Parameters
	----------
	path : str
		Path (in str) to guarantee the existence.

	Returns
	-------
	str
		The absolute path.
	"""
	if not os.path.exists(path):
		os.makedirs(path)
	return os.path.abspath(path)


def clear_dir(target: str, extension: str):
	file_pattern = f"*.{extension}"
	files_to_delete = glob.glob(os.path.join(target, file_pattern))
	for file in files_to_delete:
		os.remove(file)


def plot_dir() -> str:
	return guarantee_existence("./plots")


def velocity_dir() -> str:
	return guarantee_existence("./velocity")
