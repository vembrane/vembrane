import argparse

from .modules import filter, tablelize


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        dest="command", description="valid subcommands", required=True
    )
    filter.add_subcommmand(subparsers)
    tablelize.add_subcommmand(subparsers)

    args = parser.parse_args()
    if args.command == "filter":
        filter.execute(args)
    elif args.command == "tablelize":
        tablelize.execute(args)
    else:
        assert False
