

from config import make_logger
from config import parser_cli
from cylinder import run


logger = make_logger(
	level="INFO"
)


def main():
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

	run(lx, ly, max_t)


if __name__ == "__main__":
	main()
