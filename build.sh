#! /bin/sh

cmake -DCMAKE_BUILD_TYPE=Debug -S . -B ./build -DBUILD_TESTING=ON

cmake --build ./build
