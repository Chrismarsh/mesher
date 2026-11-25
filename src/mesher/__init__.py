"""Mesher entry point; inspection tools do not need the meshing runtime."""


# Import the CLI on invocation so inspecting a mesh does not load the GDAL, VTK and MPI meshing
# runtime just by importing the package.
def main():
    from .mesher_cli import main as mesher_main
    return mesher_main()
