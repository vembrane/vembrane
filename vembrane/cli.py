import argparse

from . import __version__
from .modules import annotate, filter, structured, table, tag


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
    filter.add_subcommmand(subparsers)
    table.add_subcommmand(subparsers)
    annotate.add_subcommmand(subparsers)
    tag.add_subcommand(subparsers)
    structured.add_subcommand(subparsers)

    args = parser.parse_args()
    if args.command == "filter":
        filter.execute(args)
    elif args.command == "table":
        table.execute(args)
    elif args.command == "annotate":
        annotate.execute(args)
    elif args.command == "tag":
        tag.execute(args)
    elif args.command == "structured":
        structured.execute(args)
    else:
        raise ValueError(f"Unknown subcommand {args.command}")
