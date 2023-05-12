#!/bin/bash


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
fi