#!/bin/bash

./gap2json.sed LayerScanningTables.gd > LayerScanningTables.json
./gap2json.sed ObliqueScanning.gd > ObliqueScanning.json
./gap2json.sed Orbits.gd > Orbits.json
./json2tex.py
