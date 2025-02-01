#!/usr/bin/env gap
# Copyright 2025 Bernard Field

Read("Scanning.gap");

orbits := List([1..48], FindLinearOrbitsFromLayerIT);

PrintTo("../TableData/Orbits.gd", orbits);
