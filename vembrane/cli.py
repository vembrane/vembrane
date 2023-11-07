import argparse

from . import __version__
from .modules import annotate, filter, table, tag


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    subparsers = parser.add_subparsers(
        dest="command",
        description="valid subcommands",
        required=True,
    )
    filter.add_subcommand(subparsers)
    table.add_subcommand(subparsers)
    annotate.add_subcommand(subparsers)
    tag.add_subcommand(subparsers)

    args = parser.parse_args()
    if args.command == "filter":
        filter.execute(args)
    elif args.command == "table":
        table.execute(args)
    elif args.command == "annotate":
        annotate.execute(args)
    elif args.command == "tag":
        tag.execute(args)
    else:
        raise ValueError(f"Unknown subcommand {args.command}")
