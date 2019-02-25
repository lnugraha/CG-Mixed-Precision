\#!/bin/bash
clear
# CLEAN the previous folder
rm -r build/

# CREATE a new folder build
mkdir build/
cd ./build

# INSIDE build folder
cmake ..; make -j16

if [[ $? -eq 0 ]]; then
	echo "Build SUCCESS"
	./HACS
else
	echo "Build FAIL";
fi;
