set -e

#setup build env
spack env create -d spack-deploy spack-deploy.yaml
spack env activate spack-deploy
spack install

#build current version
python -m build --sdist
rm -f dist/*.whl # we don't want to upload these
twine upload  dist/* --skip-existing

#clean up
spack env deactivate
rm -rf spack-deploy
rm -rf dist