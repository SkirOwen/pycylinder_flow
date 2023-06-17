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

from config import make_logger
from config import parser_cli
from cylinder import simulate


logger = make_logger(
	level="INFO"
)


def main():
	print(
		"PyCylinder_flow  Copyright (C) 2023  Owen Allemang\n"
		"This program comes with ABSOLUTELY NO WARRANTY.\n"
		"This is free software, and you are welcome to redistribute it\n"
		"under certain conditions; see `LICENSE` for details."
	)
	args = parser_cli()
	lx, ly = args.size
	max_t = args.iter
	logger.info(
		f"\n"
		f"Parameters\n"
		f"----------\n"
		f"{lx = }\n"
		f"{ly = }\n"
		f"{max_t = }\n"
	)

	simulate(lx, ly, max_t)


if __name__ == "__main__":
	main()
