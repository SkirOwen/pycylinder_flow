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