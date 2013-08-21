#!/bin/sh

g++ ccesensor2d.cpp -o ccesensor2d.exe `root-config --cflags --libs`
g++ ccesensor1d.cpp -o ccesensor1d.exe `root-config --cflags --libs`


