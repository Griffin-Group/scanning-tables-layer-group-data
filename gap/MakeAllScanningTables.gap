#!/usr/bin/env gap
# Copyright 2025 Bernard Field

# Takes about 10 minutes to run on a modern laptop.
Exec("date");

Read("Scanning.gap");

# Fix line width with SizeScreen
# So the files doesn't change every time I regenerate it.
SizeScreen([80,24]);

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

Exec("date");
