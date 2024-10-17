#!/usr/bin/env gap

# Load setting-free data
for s in ["frieze","rod","layer"] do;
  # The file is set as GAP output. This is not valid GAP input.
  file := InputTextFile(Concatenation(s,"_nosetting.gi"));
  # To make it valid GAP input, we read it in as a string
  # then turn "<data>" into  "<name> := <data>;".
  stream := InputTextString(Concatenation(s," := ",ReadAll(file),";"));
  # Then, to execute this string, we feed it into an InputTextString,
  # which is an InputStream, and Read it.
  Read(stream);
  # Make sure to close the streams when we're done.
  CloseStream(file);
  CloseStream(stream);
od;

# We now have frieze, rod, and layer in the global namespace.
# These are lists of records.
# Currently, each record has a single entry, "1", corresponding to the
# default setting.

# First, let us correct some of the default setting labels.

# All frieze groups: default setting is "a"
for i in [1..7] do;
  r := frieze[i];
  r.a := r.1;
  Unbind(r.1);
od;
# Rod groups 3 - 22: default setting is "abc"
for i in [3..22] do;
  r := rod[i];
  r.abc := r.1;
  Unbind(r.1);
od;
# Some layer groups between 8-45: default setting is "a"
ab_layer_nums := Concatenation([8..18],[20,24],[27..36],[38],[40..43],[45]);
for i in ab_layer_nums do;
  r := layer[i];
  r.a := r.1;
  Unbind(r.1);
od;
# Layer groups 52, 62, 64: default setting/origin choice is "2"
for i in [52,62,64] do;
  r := layer[i];
  r.2 := r.1;
  Unbind(r.1);
od;

# Add other settings
# Q is TransposedMat(Inverse(P)), where P is on Bilbao

ConjugateGenerators := function(r, Q, a, b)
  local gens;
  gens := List(r.(a).generators, g -> Inverse(Q) * g * Q );
  r.(b) := rec(generators := gens);
  # Optional: reduce translation to home unit cell
  # Seems to be: if translation is <0, then make it 0<x<=1.
end;

# Frieze groups, a -> b
Q := [[0,-1,0], [1,0,0], [0,0,1]];
for i in [1..7] do;
  ConjugateGenerators(frieze[i], Q, "a", "b");
od;

# Rod groups 3-22.
# abc -> bac
Q := [[0,-1,0,0], [1,0,0,0], [0,0,1,0], [0,0,0,1]];
for i in [3..22] do;
  ConjugateGenerators(rod[i], Q, "abc", "bac");
od;
# abc -> cba
Q := [[0,0,1,0], [0,1,0,0], [-1,0,0,0], [0,0,0,1]];
for i in [3..22] do;
  ConjugateGenerators(rod[i], Q, "abc", "cba");
od;
# abc -> bca
Q := [[0,0,1,0], [1,0,0,0], [0,1,0,0], [0,0,0,1]];
for i in [3..22] do;
  ConjugateGenerators(rod[i], Q, "abc", "bca");
od;
# abc -> acb
Q := [[1,0,0,0], [0,0,1,0], [0,-1,0,0], [0,0,0,1]];
for i in [3..22] do;
  ConjugateGenerators(rod[i], Q, "abc", "acb");
od;
# abc -> cab
Q := [[0,-1,0,0], [0,0,1,0], [-1,0,0,0], [0,0,0,1]];
for i in [3..22] do;
  ConjugateGenerators(rod[i], Q, "abc", "cab");
od;

# Rod group, tetragonal
Q := [[1/2,-1/2,0,0], [1/2,1/2,0,0], [0,0,1,0], [0,0,0,1]];
for i in [35,37,38,41] do;
  ConjugateGenerators(rod[i], Q, "1", "2");
od;
# Rod group, trigonal/hexagonal
Q := [[1/3,-1/3,0,0], [1/3,2/3,0,0], [0,0,1,0], [0,0,0,1]];
for i in Concatenation([46..52], [70..72], [75]) do;
  ConjugateGenerators(rod[i], Q, "1", "2");
od;

# Layer group, L5 and L7
# a -> n
Q := [[-1,-1,0,0], [1,0,0,0], [0,0,1,0], [0,0,0,1]];
for i in [5,7] do;
  ConjugateGenerators(layer[i], Q, "1", "2");
od;
# a -> b
Q := [[0,1,0,0], [-1,-1,0,0], [0,0,1,0], [0,0,0,1]];
for i in [5,7] do;
  ConjugateGenerators(layer[i], Q, "1", "3");
od;

# Layer group, orthorhombic
Q := [[0,-1,0,0], [1,0,0,0], [0,0,1,0], [0,0,0,1]];
for i in ab_layer_nums do;
  ConjugateGenerators(layer[i], Q, "a", "b");
od;

# Layer groups, different origins
Q := [[1,0,0,0], [0,1,0,0], [0,0,1,0], [1/4,1/4,0,1]];
for i in [52, 62, 64] do;
  ConjugateGenerators(layer[i], Q, "2", "1");
od;

# Pre-compute the auxiliary group information
# I won't do normalizers yet, as I don't fully grasp their consequences and uses
# But I shall do the translation basis. (This also simplifies the generators.)

SetDataBasis := function(r)
  local G, d, T;
  # Calculate the translation basis
  G := AffineCrystGroup(r.generators);
  T := TranslationBasis(G);
  # Record the translation basis
  r.basis := T;
  # Delete trivial generators
  d := Length(T[1]);
  r.generators := Filtered(r.generators, x -> x{[1..d]}{[1..d]} <> IdentityMat(d));
end;

for thelist in [frieze,rod,layer] do;
  for item in thelist do;
    for key in RecNames(item) do;
      SetDataBasis(item.(key));
    od;
  od;
od;

# Write to file
PrintTo("frieze.gi", frieze);
PrintTo("rod.gi", rod);
PrintTo("layer.gi", layer);
