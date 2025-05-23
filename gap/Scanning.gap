#!/usr/bin/env gap
# Copyright 2025 Bernard Field

LoadPackage("NumericalSgps",false); # Includes CeilingOfRational method.
#Read("RepresentativeSymmetryOps.gap");
Read("Identify_Group.gap"); # RepresentativeSymmetryOp.gap is read in this one
Read("FindOrigin.gap");

SetCrystGroupDefaultAction(LeftAction);

verbose := false;

AssertLayerGroup := function(Glayer, name)
  #% Raises error (with custom name) if Glayer is not a layer group
  if not IsAffineCrystGroup(Glayer) then;
    Error("Expected ", name, " to be an AffineCrystGroup.");
  fi;
  if Length(TranslationBasis(Glayer)) <> 2 then;
    Error("Expected ", name, " to have a 2-dimensional basis, instead got ",String(Length(TranslationBasis(Glayer))));
  fi;
  if Length(InternalBasis(Glayer)) <> 3 then;
    Error("Expected ", name, " to act in 3-dimensional space, instead got ",String(Length(InternalBasis(Glayer))));
  fi;
end;


# Can I throw all of these together into some self-contained method?
SectionalGroupAndTransformationOfLayerGroupNoOrigin := function(Glayer, vrod, vscan, s)
  local vz, vz2, vscan2, vrod2, transformation, Gscan, GscanOps, GrodOps, g, i, Grod, shiftMat;
  #% Given a layer group, direction for the rod, direction for scanning, and 
  #% position s along scanning vector, return a rod group (in basis [z, scan, rod])
  #% and the transformation matrix from the old to new basis.
  if verbose then; Print("SectionalGroupAndTransformationOfLayerGroupNoOrigin called.\n"); fi;
  # Verify arguments
  AssertLayerGroup(Glayer, "Glayer");
  if not IsVector(vrod) then;
    Error("Expected vrod to be a vector.");
  fi;
  if not IsVector(vscan) then;
    Error("Expected vscan to be a vector.");
  fi;
  if not IsRat(s) then;
    Error("Expected s to be a rational number.");
  fi;
  # First step...
  # Find the out-of-plane element
  vz := Filtered(SymmetricInternalBasis(Glayer), v -> not v in TranslationBasis(Glayer))[1];
  # Ensure a right-handed coordinate system.
  if Determinant([vscan, vz, vrod]) < 0 then
    vz := -vz;
  fi;
  if verbose then; Print("vz = ",String(vz),"\n"); fi;
  # Next step, transform to the scanning basis, with origin shift
  # Let the order by vz, vscan, vrod
  vz2 := [0,1,0];
  vscan2 := [1,0,0];
  vrod2 := [0,0,1];
  transformation := TranslationAffineMatrix( s * vscan2 );
  transformation{[1..3]}{[1..3]} := Inverse(TransposedMat([vscan, vz, vrod]));
  if verbose then; Print("Transformation matrix:\n"); Display(transformation); fi;
  Gscan := \^(Glayer, transformation);
  # We could call findRodGroupFromScanning, but that assumes we've filtered out
  # the point group, and haven't done the origin shift. So I'll do it myself.
  GscanOps := RepresentativeSymmetryOps(Gscan, [vrod2, vscan2]);
  if verbose then;
    Print("Representative operations of transformed layer group:\n");
    for g in GscanOps do; Display(g); od;
  fi;
  # Filter the operations to those with appropriate translation and rotation
  # Integer along scanning direction
  GrodOps := Filtered(GscanOps, g -> IsInt(g[1][4]) and
	     (g{[1..3]}{[1..3]} * vrod2 = vrod2 or g{[1..3]}{[1..3]} * vrod2 = -vrod2));
  # Enforce 0 mod 1 -> 0.
  for i in [1..Length(GrodOps)] do;
    g := MutableCopyMatrix(GrodOps[i]);
    g[1][4] := 0;
    GrodOps[i] := g;
  od;
  if verbose then;
    Print("Representative operations of rod group:\n");
    for g in GrodOps do; Display(g); od;
  fi;
  # Add the lattice translation.
  Add(GrodOps, TranslationAffineMatrix(vrod2));
  # Create the group
  Grod := AffineCrystGroup(GrodOps);
  return [Grod, transformation];
end;

SectionalGroupAndTransformationOfLayerGroup := function(Glayer, vrod, vscan, s)
  local res, Grod, transformation, shiftMat;
  #% Given a layer group, direction for the rod, direction for scanning, and 
  #% position s along scanning vector, return a rod group (in basis [z, scan, rod])
  #% and the transformation matrix from the old to new basis. Also transforms to standard origin.
  res := SectionalGroupAndTransformationOfLayerGroupNoOrigin(Glayer, vrod, vscan, s);
  Grod := res[1];
  transformation := res[2];
  # Shift to standard origin
  shiftMat := StandardOriginTransformationOfGroup(Grod);
  Grod := Grod ^ shiftMat;
  transformation := shiftMat * transformation;
  if verbose then;
    Print("Translated origin by ",String(shiftMat{[1..3]}[4]),"\n");
  fi;
  return [Grod, transformation];
end;

SectionalGroupOfLayerGroup := function(Glayer, vrod, vscan, s)
  #% Given a layer group, direction for the rod, direction for scanning, and 
  #% position s along scanning vector, return a rod group (in basis [scan, z, rod]).
  return SectionalGroupAndTransformationOfLayerGroup(Glayer, vrod, vscan, s)[1];
end;

FloorOfRational := r -> -CeilingOfRational(-r);

RoundRational := function(r)
  local fr, cr;
  if IsInt(r) then
    return r;
  else
    fr := FloorOfRational(r);
    cr := CeilingOfRational(r);
    if r - fr > cr - r then
      return cr;
    else
      return fr;
    fi;
  fi;
end;

ScanningVectorFromIntVector := function(vec)
  local p0, q0, sol, nullvec, x;
  #% Helper function for ScanningTable. Gets a scanning vector from the rod vector vec, which must be an integer 3-vector with zero third component.
  if not ForAll(vec, IsInt) then
    Error("vec is not a crystallogrpahic vector in standard basis! ", vec);
  fi;
  # If vec = [m, n, 0] and output is vscan = [p, q, 0], then we solve
  # mq - np = 1 (i.e. vec x vscan = 1). This ensures our unit cell is
  # of unit area and the coordinate system is right-handed.
  sol := SolutionNullspaceIntMat([[vec[1]],[-vec[2]]], [1]);
  p0 := sol[1][2];
  q0 := sol[1][1];
  nullvec := [sol[2][1][2], sol[2][1][1], 0];
  # We then find the solution closest to being right-handed.
  # x is how much nullvec to add to the specific solution [p0,q0,0] to
  # minimise the dot product squared between vec and vscan.
  x := -vec * [p0,q0,0] / (nullvec * vec);
  # Round to the nearest integer.
  x := RoundRational(x);
  return [p0, q0, 0] + nullvec * x;
end;

ScanningVectorFromIntVectorAndGroup := function(vec, G)
  local OD, gens, daux, vaux;
  # Use group information to choose the ideal scanning vector giving a conventional basis
  OD := List(GeneratorsOfGroup(PointGroup(G)), Order);
  if ForAny(OD, x -> x >= 3) then
    # Trigonal, tetragonal, or hexagonal.
    # Catch cases where there are no in-plane operations.
    # In these cases, the scanning group is not centering, so we should pick a primitive
    # basis as our conventional basis.
    gens := Filtered(GeneratorsOfGroup(G),
            x -> x{[1..3]}{[1..3]} <> IdentityMat(3) and x{[1..3]}{[1..3]} <> - IdentityMat(3));
    if ForAll(List(gens, AxisOfOperator), x -> x = [0,0,1]) then
      return ScanningVectorFromIntVector(vec);
    fi;
  fi;
  # Tetragonal
  if 4 in OD then
    if vec = [1,1,0] then
      return [-1,1,0];
    elif vec = [1,-1,0] then
      return [1,1,0];
    else
      return ScanningVectorFromIntVector(vec);
    fi;
  # Trigonal/hexagonal
  elif 3 in OD or 6 in OD then
    if vec = [1,0,0] then
      return [1,2,0];
    elif vec = [0,1,0] then
      return [-2,-1,0];
    elif vec = [1,1,0] then
      return [-1,1,0];
    elif vec = [1,2,0] then
      return [-1,0,0];
    elif vec = [2,1,0] then
      return [0,1,0];
    elif vec = [1,-1,0] then
      return [1,1,0];
    else
      return ScanningVectorFromIntVector(vec);
    fi;
  else
    # Catch the case of inclined scanning in centering groups.
    if ForAny(vec, x -> AbsInt(x) > 1) and
      (TranslationAffineMatrix([1/2,1/2,0]) in G) then
      # Solve an auxiliary problem
      vaux := [-vec[1] - vec[2], vec[2] - vec[1], 0];
      daux := ScanningVectorFromIntVector(vaux);
      # Get the scanning vector to be a centering vector
      return [(daux[1]+daux[2])/2, (daux[1]-daux[2])/2, 0];
    else
      return ScanningVectorFromIntVector(vec);
    fi;
  fi;
end;

# Get the scanning group (I don't use it for my calculations, but it's reported in the tables)
# This is the equitranslational subgroup of G whose point group contains all the elements of G
# which leave the plane/rod (defined by vrod) invariant.
# By the scanning theorem, the scanning table of G is identical to that of its scanning group.
ScanningGroupOfLayerGroup := function(G, args...)
  local vrod, vscan, action, ops, g, M, t, H, vz, P;
  #% Get the scanning group of a layer group
  # Unpack the arguments
  vrod := args[1];
  if Length(args) >= 2 then
    vscan := args[2];
  else
    # Get default scanning vector (for getting the right setting).
    vscan := ScanningVectorFromIntVectorAndGroup(vrod,G);
  fi;
  # Determine if G is left or right action
  if IsAffineCrystGroupOnLeft(G) then
    if IsAffineCrystGroupOnRight(G) then
      action := CrystGroupDefaultAction;
    else
      action := LeftAction;
    fi;
  elif IsAffineCrystGroupOnRight(G) then
    action := RightAction;
  else
    ErrorNoReturn("G is not an AffineCrystGroup!\n");
  fi;
  if action = LeftAction then
    G := TransposedMatrixGroup(G);
  fi;
  # Eliminate the elements of G whose linear part changes +-vrod.
  ops := [];
  for g in RepresentativeSymmetryOps(G) do
    M := g{[1..3]}{[1..3]};
    if vrod * M = vrod or vrod * M = -vrod then
      Add(ops, g);
    fi;
  od;
  # Add the translations of G back into the mix
  for t in TranslationBasis(G) do
    g := IdentityMat(4);
    g[4]{[1..3]} := t;
    Add(ops, g);
  od;
  # Create the scanning group
  H := AffineCrystGroupOnRight(ops);
  # Transform into the scanning basis
  vz := Filtered(SymmetricInternalBasis(G), v -> not v in TranslationBasis(G))[1];
  P := IdentityMat(4);
  P{[1..3]}{[1..3]} := [vrod, vscan, vz];
  H := H ^ Inverse(P);
  if action = LeftAction then
    H := TransposedMatrixGroup(H);
  fi;
  return H;
end;

# Now that I can do that, let's find all special values of s.
# A time-efficient implementation might do some clever caching to return the
# groups with the values of s. But it would be far simpler to just find all
# the special values of s then create groups using the above method.

SpecialScanSOfLayerGroup := function(Glayer, vrod, vscan)
  local vz, transformation, vscan2, vrod2, vz2, Gscan, g, gR, q, t, n, sSet;
  #% For a layer group, center line direction, and scanning direction,
  #% returns a Set of numbers between 0 and 1 where special scanning positions
  #% are. vrod and vscan must be crystallographic vectors (i.e. integer
  #% components).
  # Validate arguments
  AssertLayerGroup(Glayer, "Glayer");
  if not IsVector(vrod) then;
    Error("Expected vrod to be a vector.");
  fi;
  if not IsVector(vscan) then;
    Error("Expected vscan to be a vector.");
  fi;
  # Find the out-of-plane element
  vz := Filtered(SymmetricInternalBasis(Glayer), v -> not v in TranslationBasis(Glayer))[1]; # Hopefully this is robust...
  # Next step, transform to the scanning basis, with origin shift
  # Let the order by vz, vscan, vrod
  vz2 := [0,1,0];
  vscan2 := [1,0,0];
  vrod2 := [0,0,1];
  transformation := IdentityMat(4);
  transformation{[1..3]}{[1..3]} := Inverse(TransposedMat([vscan, vz, vrod]));
  Gscan := \^(Glayer, transformation);
  # Now, scanning.
  # The center line shifts by s * scand, where s is a real number.
  # The Origin shifts by this as well, so the line is on the origin.
  # For each operation, we want to find the values of s for which the translation
  # component (t + s*(d-R*d))[idxscan] = 0 mod 1, where s*(d-R*d) is the origin shift, t is the translation component, and R is the linear component.
  # If (d-R*d) = 0, then this is only true if t = 0 mod 1, and the operation
  # belongs to the general case.
  # Otherwise, q:=(d-R*d)[idxscan], we get a special case, where for some discrete
  # values of s the operation is allowed.
  # s = n/q - t[idxscan]/q, where n is an integer.
  # We restrict ourselves to values of s between 0 and 1.
  # Once we have assigned values of valid s to each symmetry operation,
  # including for the translation subgroup, then we can build a new group
  # for each case of s.
  # So, go over each operation, extract the linear (R) and translation (t) parts,
  # get the allowed s = n/q - t/q. Make a set of special s values.
  sSet := Set([]);
  for g in RepresentativeSymmetryOps(Gscan, [vscan2, vrod2]) do;
    # We need to exclude any operations for which the point group is wrong
    gR := g{[1..3]}{[1..3]}; # Linear part of g
    if (gR * vrod2 = vrod2 or gR * vrod2 = -vrod2) then;
      q := (vscan2 - gR * vscan2)[1]; # Unit origin shift in scan coordinate
      if q <> 0 then;
        # Only have special s when q <> 0.
        t := g[1][4]; # Scanning component of translation.
        # We want to restrict s to be 0 <= s < 1.
        # n are integers selected to maintain this condition.
        if q > 0 then;
          n := [CeilingOfRational(t) .. (CeilingOfRational(q+t)-1)];
        else;
          n := [(FloorOfRational(q+t)+1) .. FloorOfRational(t)];
        fi;
        UniteSet(sSet, (n-t)/q);
      fi;
    fi;
  od;
  return sSet;
end;

# Next task: create a scanning table.

ScanningTableFunAllInfo := function(G, vrod, vscan)
  local sSet, genS, s, num, r, rS, found, len, res, setting, origin;
  # Get special s.
  sSet := SpecialScanSOfLayerGroup(G, vrod, vscan);
  # Get a general s
  if Length(sSet) = 0 then
    genS := 0;
  elif Length(sSet) = 1 then
    if 0 in sSet then
      genS := 1/2;
    else
      genS := 0;
    fi;
  else
    # Choose genS to be somewhere between two specials.
    # Make it a non-even number so it's more obvious if a property is s-dependent.
    genS := (3*sSet[1] + 4*sSet[2])/7;
  fi;
  if genS in sSet then
    Error("Somehow the general S", genS, "is in the special s set", sSet);
  fi;
  # Make a record with the IT number of the groups
  r := rec();
  r.d := vscan; # Scanning vector
  # General position: rod group number
  res := SectionalGroupAndTransformationOfLayerGroupNoOrigin(G, vrod, vscan, genS);
  r.general := ITNumberOfRodGroup(res[1]);
  r.general_group := res[1];
  r.general_setting := SettingOfRodGroupNC(res[1], r.general);
  r.general_transform := res[2];
  r.general_length := VectorLengthRatio(vrod, TranslationBasis(res[1]^Inverse(res[2]))[1]);
  r.general_s := genS;
  r.general_origin := OriginLocationOfGroup(r.general_group);
  # Special positions
  r.special := [];
  for s in sSet do
    res := SectionalGroupAndTransformationOfLayerGroupNoOrigin(G, vrod, vscan, s);
    num := ITNumberOfRodGroup(res[1]);
    len := VectorLengthRatio(vrod, TranslationBasis(res[1]^Inverse(res[2]))[1]);
    setting := SettingOfRodGroupNC(res[1], num);
    origin := OriginLocationOfGroup(res[1]);
    # Check if already present
    found := false;
    for rS in r.special do
      if rS.number = num and rS.length = len and rS.setting = setting and rS.origin = origin then
        Add(rS.s, s);
        Add(rS.transform, res[2]);
        Add(rS.group, res[1]);
        found := true;
        break;
      fi;
    od;
    if not found then
      Add(r.special, rec(number := num, length := len, s := [s],
                         group := [res[1]], transform := [res[2]],
                                  setting := setting, origin := origin)
         );
      # 'transform' varies by a translation s.
      # 'group' is often the same, but not always. It might have something to do with orbits? Or maybe setting.
    fi;
  od;
  # Get the scanning group
  r.H_group := ScanningGroupOfLayerGroup(G, vrod, vscan);
  r.H := ITNumberOfLayerGroup(r.H_group);
  r.H_setting := SettingOfLayerGroupNC(r.H_group, r.H);
  r.H_origin := OriginLocationOfGroup(r.H_group);
  return r;
end;

ScanningTableFun := function(G, vrod, vscan)
  local r, rS;
  # Delete the extraneous group and transformation information that we don't
  # need for a simple scanning table. (It exists 
  r := ScanningTableFunAllInfo(G, vrod, vscan);
  Unbind(r.general_group);
  Unbind(r.general_transform);
  Unbind(r.general_s);
  for rS in r.special do
    Unbind(rS.group);
    Unbind(rS.transform);
  od;
  Unbind(r.H_group);
  return r;
end;

ScanningTable := function(G, vrod)
  local T, vscan, v, r;
  AssertLayerGroup(G, "G");
  # Move to standard basis
  T := Inverse(TransposedMat(SymmetricInternalBasis(G)));
  G := G ^ AugmentedMatrix(T, [0,0,0]);
  vrod := T * vrod;
  # Create vscan, a vector in TranslationBasis(G) not parallel to vrod.
  vscan := ScanningVectorFromIntVectorAndGroup(vrod,G);
  # Get the scanning table
  r := ScanningTableFun(G, vrod, vscan);
  # Convert the scanning vector back to the original basis
  r.d := Inverse(T) * r.d;
  return r;
end;

ScanningTableFromLayerIT := function(num, vrod)
  local G, vscan;
  # Because IT gives it in conventional basis, I don't need to transform to
  # standard form.
  G := LayerGroupIT(num);
  vscan := ScanningVectorFromIntVectorAndGroup(vrod,G);
  return ScanningTableFun(G, vrod, vscan);
end;

GetRodDirections := function(num)
  local rods;
  # Iterate over key rod vectors for a IT Number
  rods := [[1,0,0]];
  if num >= 8 or num = 5 or num = 7 then
    # Beyond oblique, or p11a or p112/a
    Add(rods, [0,1,0]);
  fi;
  if num = 7 or num = 5 then
    # A very special position for p112/a (or, rather, p112/n)
    # and p11a (or p11n)
    Add(rods, fail);
  fi;
  if num >= 49 then
    # Square or hexagonal
    Add(rods, [1,1,0]);
    Add(rods, [1,-1,0]);
  fi;
  if num >= 65 then
    # Hexagonal
    Add(rods, [1,2,0]);
    Add(rods, [2,1,0]);
  fi;
  return rods;
end;

FullScanningTableFromLayerIT := function(num)
  local table, r, vrod, found, x;
  table := [];
  # Iterate over all the distinct rod vectors up to some distance
  for vrod in GetRodDirections(num) do
    Print("Scanning layer # ", num, " and orientation ", vrod, ".\n");
    if vrod = fail then
      # Very special case for p112/n (L7, setting "2").
      vrod := [0,1,0];
      r := ScanningTableFun(LayerGroupIT(num), vrod, [-1,-1,0]);
    else
      r := ScanningTableFromLayerIT(num, vrod);
    fi;
    r.orientation := [vrod];
    r.d := [r.d];
    # Consolidate, such that each distinct set of rod groups appears only once,
    # and we report all the orientations that give the same result together
    found := false;
    for x in table do
      # Must do IsRat comparison before value equality because, in some cases,
      # general_length may be a float (because it's irrational).
      if x.general = r.general and
        IsRat(x.general_length) = IsRat(r.general_length) and
        x.general_length = r.general_length and
        x.general_setting = r.general_setting and
        x.general_origin = r.general_origin and
        x.special = r.special and
        x.H = r.H and
        x.H_setting = r.H_setting and
        x.H_origin = r.H_origin then
        found := true;
        Append(x.d, r.d);
        Append(x.orientation, r.orientation);
        break;
      fi;
    od;
    if not found then
      Add(table, r);
    fi;
  od;
  return table;
end;

# As a supplemental function (which we shall perform separately from the rest
# of the tables for now), find the linear orbits of each of the scanning groups.

FindLinearOrbitsNC := function(G, vrod, vscan, origin)
  local d, Ss, orbits, ops, proj, s1, s2, P1, P2, myorbit, g;
  # N.B. Algorithm will FAIL for layer groups greater than 48.
  # Convert to augmented vectors if necessary
  d := DimensionOfMatrixGroup(G) - 1;
  if Length(vrod) = d then
    vrod := Concatenation(vrod, [0]);
  fi;
  if Length(vscan) = d then
    vscan := Concatenation(vscan, [0]);
  fi;
  if Length(origin) = d then
    origin := Concatenation(origin, [1]);
  fi;
  # Get the Ss we want to check.
  Ss := [0,1/4,1/2,3/4];
  # Initialise list of orbits.
  orbits := [];
  # Get operators.
  # Specify the conventional translation basis
  ops := RepresentativeSymmetryOps(G, [[1,0,0],[0,1,0]]);
  # Get to consistent action.
  if IsAffineCrystGroupOnRight(G) then
    ops := List(ops, TransposedMat);
  fi;
  # If we're in a centering group, add the centering operation
  #if [1/2,1/2,0] in TranslationBasis(G) then
  #  Add(ops, [[1,0,0,1/2],[0,1,0,1/2],[0,0,1,0],[0,0,0,1]]);
  #fi;
  # I've done some maths. There is a vector we dot product by to project
  # onto the original line.
  proj := ( (vrod*vrod) * vscan - (vrod*vscan) * vrod ) / ( (vrod*vrod)*(vscan*vscan) - (vrod*vscan)^2 );
  # Iterate over each s.
  for s1 in Ss do
    myorbit := [s1];
    # Get the starting point.
    P1 := origin + s1 * vscan;
    # Iterate over each operator
    for g in ops do
      P2 := g * P1;
      s2 := (P2 - origin) * proj;
      # Round
      s2 := s2 - FloorOfRational(s2);
      if s2 <> s1 and s2 in Ss then
        # Pop s2 from Ss so we don't duplicate.
        # (But don't remove current element, because then the for loop will get confused.)
        Remove(Ss, Position(Ss, s2));
      fi;
      if not s2 in myorbit then
        Add(myorbit, s2);
      fi;
    od;
    Add(orbits, myorbit);
  od;
  # Get the general orbit
  s1 := 1/7;
  myorbit := [s1];
  # Get the starting point.
  P1 := origin + s1 * vscan;
  # Iterate over each operator
  for g in ops do
    P2 := g * P1;
    s2 := (P2 - origin) * proj;
    # Round
    s2 := s2 - FloorOfRational(s2);
    if not s2 in myorbit then
      Add(myorbit, s2);
    fi;
  od;
  # Sort, for consistency
  Add(orbits, Set(myorbit));
  return orbits;
end;

FindLinearOrbitsFromLayerIT := function(num)
  local G, r, vrod, vscan, Ss, orbits;
  if num > 48 then
    Error("FindLinearOrbits only works for Layer groups up to 48.");
  fi;
  # Set up things
  G := LayerGroupIT(num);
  r := rec();
  r.a := FindLinearOrbitsNC(G, [1,0,0], [0,1,0], [0,0,0]);
  r.b := FindLinearOrbitsNC(G, [0,1,0], [-1,0,0], [0,0,0]);
  if [1/2,1/2,0] in TranslationBasis(G) then
    # Centering group, so origin shifts should also be checked.
    r.ao := FindLinearOrbitsNC(G, [1,0,0], [0,1,0], [1/4,1/4,0]);
    r.bo := FindLinearOrbitsNC(G, [0,1,0], [-1,0,0], [1/4,1/4,0]);
  fi;
  return r;
end;
