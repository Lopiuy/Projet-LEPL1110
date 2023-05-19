#!/bin/bash

cd ..

if [ "$#" == "0" ]; then
  cd ProjectPreProcessor/build/ || echo "Error: ProjectPreProcessor/build/ does not exist"
  make
  ./myFem maison
  cd ../../..
  cd build/ || echo "Error: Project/build/ does not exist"
  make
  ./myFem
  cd ../src
  cd ProjectPostProcessor/build/ || echo "Error: ProjectPostProcessor/build/ does not exist"
  make
  ./myFem
  cd ../..

elif [ "$1" == "pre" ]; then
  cd ProjectPreProcessor/build/ || echo "Error: ProjectPreProcessor/build/ does not exist"
  make
  ./myFem maison view
  cd ../..

fi
