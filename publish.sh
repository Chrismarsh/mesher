set -e
# spack env create . py-mesher-spack.yaml
spack env activate . --sh
# spack install
python -m build --sdist
rm -f dist/*.whl # we don't want to upload these
twine upload  dist/* --skip-existing
