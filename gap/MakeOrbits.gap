#!/usr/bin/env gap

Read("Scanning.gap");

orbits := List([1..48], FindLinearOrbitsFromLayerIT);

PrintTo("../TableData/Orbits.gd", orbits);
