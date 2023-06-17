# Copyright (C) 2023  Owen Allemang

# This file is part of PyCylinder_flow.
# PyCylinder_flow is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.

# PyCylinder_flow is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with Foobar.
# If not, see <https://www.gnu.org/licenses/>.

import logging

import argparse
from argparse import Namespace


def make_logger(level: str = "INFO") -> logging.Logger:
	FORMAT = "%(message)s"
	logging.basicConfig(
		level=logging.WARNING, format=FORMAT, datefmt="[%X]"
	)

	logger = logging.getLogger("rdp")
	logger.setLevel(level)
	return logger


def parser_cli() -> Namespace:
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"-i", "--iter",
		type=int,
		help="Number of iteration to simulate"
	)

	parser.add_argument(
		"-s", "--size",
		nargs=2,
		type=int,
		help="Size x, and y of the simulation"
	)

	args = parser.parse_args()
	return args