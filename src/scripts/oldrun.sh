#!/bin/bash

cd ..

if [ "$#" == "0" ]; then
  cd ProjectPreProcessor/build/
  make
  ./myFem
  cd ../..
  cd Project/build/
  make
  ./myFem
  cd ../..
  cd ProjectPostProcessor/build/
  make
  ./myFem
  cd ../..

elif [ "$1" == "pre" ]; then
  cd ProjectPreProcessor/build/
  make
  ./myFem view
  cd ../..

elif [ "$1" == "proj" ]; then
  cd Project/build/
  make
  ./myFem
  cd ../..


elif [ "$1" == "post" ]; then
  cd ProjectPostProcessor/build/
  make
  ./myFem
  cd ../..

elif [ "$1" == "projdb" ]; then
  cd Project/build/
  make
  gdb ./myFem
  cd ../..

elif [ "$1" == "predb" ]; then
  cd ProjectPreProcessor/build/
  make
  gdb ./myFem
  cd ../..

fi