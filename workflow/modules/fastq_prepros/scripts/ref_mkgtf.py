import argparse
import subprocess


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", required=True)
    p.add_argument("--filtered-gtf", required=True)
    args = p.parse_args()

    subprocess.run(["celescope", "utils", "mkgtf", args.gtf, args.filtered_gtf], check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

