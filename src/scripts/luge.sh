#!/bin/bash

cd ..

if [ "$#" == "0" ]; then
  cd ProjectPreProcessor/build/ || echo "Error: ProjectPreProcessor/build/ does not exist"
  make
  ./myFem luge
  cd ../../..
  cd build/ || echo "Error: Project/build/ does not exist"
  make
  ./myFem luge
  cd ../src
  cd ProjectPostProcessor/build/ || echo "Error: ProjectPostProcessor/build/ does not exist"
  make
  ./myFem
  cd ../..

elif [ "$1" == "pre" ]; then
  cd ProjectPreProcessor/build/ || echo "Error: ProjectPreProcessor/build/ does not exist"
  make
  ./myFem luge view
  cd ../..

fi
