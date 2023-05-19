#!/bin/bash

cd scripts/ || echo "Error: scripts/ does not exist"

if [ "$#" == "0" ]; then
  ./luge.sh

elif [ "$1" == "luge" ]; then
  if [ "$2" == "form" ]; then
    ./luge.sh pre
  else
    ./luge.sh
  fi

elif [ "$1" == "basic" ]; then
  if [ "$2" == "form" ]; then
    ./basic.sh pre
  else
    ./basic.sh
  fi

elif [ "$1" == "poutre" ]; then
  if [ "$2" == "form" ]; then
    ./poutre.sh pre
  else
    ./poutre.sh
  fi


elif [ "$1" == "maison" ]; then
  if [ "$2" == "form" ]; then
    ./maison.sh pre
  else
    ./maison.sh
  fi

elif [ "$1" == "build" ]; then
  ./build.sh

fi
