#!/bin/sh

g++ ccesensor2d.C -o ccesensor2d.exe `root-config --cflags --libs`
g++ ccesensor1d.C -o ccesensor1d.exe `root-config --cflags --libs`


