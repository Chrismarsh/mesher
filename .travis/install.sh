#!/bin/bash

set -e
set -x

if [ "$TRAVIS_OS_NAME" = "osx" ]; then

    bash -c 'echo $pyv'
    bash -c 'echo $TRAVIS_PYTHON_VERSION'
    bash -c 'echo $TRAVIS_OS_NAME'

    brew update || true
    brew upgrade || true
    brew install gdal || brew upgrade gdal
    brew outdated pyenv || brew upgrade pyenv
    brew install pyenv-virtualenv

    eval "$(pyenv init -)"
    pyenv install $pyv
    pyenv virtualenv $pyv mesher
    pyenv rehash
    pyenv activate mesher

else
  sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y;
  sudo apt-get update -qq
  sudo apt-get install g++-7
  sudo apt-get install libgdal-dev
  sudo apt-get install python-gdal
  sudo apt-get install gdal-bin
#  eval "CC=gcc-7 && CXX=g++-7"
fi
