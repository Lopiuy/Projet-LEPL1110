#!/bin/bash

cd ..

cd ProjectPreProcessor/build/ || echo "Error: ProjectPreProcessor/build/ does not exist"
cmake ..
cd ../../..

cd build/ || echo "Error: Project/build/ does not exist"
cmake ..
cd ../src

cd ProjectPostProcessor/build/ || echo "Error: ProjectPostProcessor/build/ does not exist"
cmake ..
cd ../..
