import argparse
from pathlib import Path


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output", required=True)
    args = p.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("MT-\nmt-\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

