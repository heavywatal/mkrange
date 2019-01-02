#!/bin/sh
for x in {6000..6007}; do
  echo "\n## inputfile$x\n"
  mkdir -p output$x
  cd output$x
  ln -sf ../inputfile$x Inputfile
  ../../a.out steep
  cd ..
done
