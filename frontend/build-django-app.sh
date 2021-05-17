#!/bin/bash

mkdir -p gibbsml
cp -pr ../LICENSE ../MANIFEST.in ../README.md \
       ../requirements.txt ../setup.py ../gibbsml \
       ./gibbsml/

sudo docker build -t batterycycling .
