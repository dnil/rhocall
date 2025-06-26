import sys

from rhocall.cli import cli

if __name__ == "__main__":
    # exit using whatever exit code the CLI returned
    sys.exit(cli())
