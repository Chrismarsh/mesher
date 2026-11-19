"""Mesher entry point; inspection tools do not need the meshing runtime."""


def main():
    from .mesher_cli import main as mesher_main
    return mesher_main()
