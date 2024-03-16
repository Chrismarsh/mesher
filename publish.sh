python -m build
rm -f dist/*.whl # we don't want to upload these
twine upload  dist/*
