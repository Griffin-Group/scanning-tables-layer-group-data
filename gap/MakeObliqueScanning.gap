#!/usr/bin/env gap

# Find the oblique scanning group of each layer group.

Read("Scanning.gap");

groups := List([1..80], LayerGroupIT);
scanning := List(groups, g -> ITNumberOfLayerGroup(ScanningGroupOfLayerGroup(g, [3,-2,0])));

PrintTo("../TableData/ObliqueScanning.gd", scanning);
