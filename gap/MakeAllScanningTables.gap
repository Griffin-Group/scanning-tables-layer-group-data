#!/usr/bin/env gap
# Copyright 2024 Bernard Field

Read("Scanning.gap");

fname := "../TableData/LayerScanningTables.gd";
PrintTo(fname, "tables := [\n");

tables := [];;

for i in [1..80] do
  table := FullScanningTableFromLayerIT(i);;
  Add(tables, table);
  AppendTo(fname, table);
  if i < 80 then
    AppendTo(fname, ",\n");
  else
    AppendTo(fname, "\n];\n");
  fi;
od;
