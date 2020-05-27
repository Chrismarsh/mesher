#!/bin/bash

set -e

if [ "$TRAVIS_OS_NAME" = "osx" ]; then
    eval "$(pyenv init -)"
    pyenv activate mesher
fi

if [ "$test_conda" = "1" ]; then exit(0) ; fi

pip install twine
pip install conan
pip install scikit-build==0.10.0
pip install ninja
pip install wheel

if [ "$TRAVIS_OS_NAME" = "osx" ]; then
  python setup.py sdist bdist_wheel
else
  python setup.py sdist #no binary wheels on linux at the moment as we link against non PEP0513 .so
  #https://www.python.org/dev/peps/pep-0513/
fi
twine upload  --skip-existing dist/*