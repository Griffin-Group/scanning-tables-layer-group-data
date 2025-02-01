#!/bin/env gap
# Copyright 2025 Bernard Field

LoadPackage("NumericalSgps", false); # For CeilingOfRational
Read("RepresentativeSymmetryOps.gap");

# Here we want to take a random rod group and identify its International Tables number

# Step 1: determine the point group.
# This can be done in a basis-independent way by considering the determinant
# and order of each operator, then comparing.
# Ah, but also need to check for inversion separately.

OrderDetOfMatrix := function(M)
  #% To distinguish from mirror, inversion has order -1.
  if M = -Identity(M) then
    return [-1, -1];
  else
    return [Order(M), Determinant(M)];
  fi;
end;

OrderDetSymOfMatrix := function(M)
  local d, OD, P, eigs;
  if not IsMatrix(M) then
    Error("M is not a matrix.");
  fi;
  d := Length(M) - 1; # Dimension
  P := M{[1..d]}{[1..d]}; # Linear part of matrix
  OD := OrderDetOfMatrix(P);
  # An operation is symmorphic if it has a fixed point,
  # that is, an eigenvector which represents a point [x,y,z,1].
  # (You could also do it by counting linearly independent eigenvectors,
  # but only if you can compute complex eigenvectors, which GAP can't.)
  if IsAffineMatrixOnLeft(M) then
    M := TransposedMat(M); # Eigenvectors acts on the right.
  fi;
  eigs := Eigenvectors(Rationals, M);
  Add(OD, ForAny(eigs, x -> x[d+1] <> 0));
  return OD;
end;
OrderDetSymMatrix := function(M)
  local res;
  res := OrderDetSymOfMatrix(M);
  Add(res, M);
  return res;
end;

DecomposeMatrixTranslationOnRight := function(M)
  local d, W, w, n, wg, wl;
  #% Input: AffineMatrixOnRight. Out: intrinsic translation, location translation
  if not IsAffineMatrixOnRight(M) then
    ErrorNoReturn("M is not an AffineMatrixOnRight");
  fi;
  # Dimension
  d := Length(M) - 1;
  # Linear component
  W := M{[1..d]}{[1..d]};
  # Translation component
  w := M[d+1]{[1..d]};
  # Order
  n := Order(W);
  # Intrinsic part of the translation can be found by raising M to power n
  # The location part doesn't accumulate over repeated applications (and 
  # vanishes for the identity), while the intrinsic part remains and accumulates
  wg := (M^n)[d+1]{[1..d]} / n;
  # The location part is found by simple subtraction
  wl := w - wg;
  # Return
  return [wg, wl];
end;
IsMatrixSymmorphicEquivalent := function(M, G)
  local T, wg, Tg, d, W, t, v, L, dividezeronan, i, ns;
  #% Operation M in G. If we translate it, does it become symmorphic?
  # N.B. we only consider the primitive cell, not the conventional cell.
  if not M in G then
    Error("Matrix is not in G");
  fi;
  # Go to RightAction
  if IsAffineCrystGroupOnLeft(G) then
    M := TransposedMat(M);
  fi;
  # Translation basis
  T := TranslationBasis(G);
  # Get intrinsic part of translation
  wg := DecomposeMatrixTranslationOnRight(M)[1];
  # Catch trivial case
  if wg = Zero(wg) then
    return true;
  fi;
  if Length(T) = 0 then
    # wg is non-zero if we got here.
    return false;
  fi;
  # Now for the hard part. We need to determine whether a lattice translation,
  # applied to M, can make wg zero. (Not an origin shift, but an extra 
  # operation.)
  # Get how the translation basis will affect wg.
  Tg := [];
  d := Length(M)-1; # Dimension
  W := M{[1..d]}{[1..d]}; # Linear part
  for t in T do
    Add(Tg, DecomposeMatrixTranslationOnRight(AugmentedMatrix(W, t))[1]);
  od;
  # Now we solve the equation wg = TransposedMat(Tg) * x for x.
  # Or, at least, we test whether x has integer solutions.
  return IntSolutionMat(Tg, wg) <> fail;
end;

OrderDetMatrixOfGroup := function(G)
  local ops, OD, i;
  ops := RepresentativeSymmetryOps(G);
  # OrderDetOfMatrix on linear part of matrix
  OD := List(ops, x -> OrderDetOfMatrix(x{[1..Length(x)-1]}{[1..Length(x)-1]}));
  for i in [1..Length(ops)] do
    Add(OD[i], ops[i]);
  od;
  return OD;
end;

OrderDetSymMatrixOfGroup := function(G)
  local ops, ODS, elem;
  ops := RepresentativeSymmetryOps(G);
  ODS := List(ops, OrderDetSymMatrix);
  # We now have a special task. Sometimes, RepresentativeSymmetryOps might
  # return a matrix which isn't symmorphic when another representative is.
  # We need to catch this case.
  for elem in ODS do
    if not elem[3] then
      # Reported as not symmorphic
      elem[3] := IsMatrixSymmorphicEquivalent(elem[4], G);
    fi;
  od;
  Sort(ODS);
  return ODS;
end;
OrderDetSymOfGroup := function(G)
  local ODS;
  ODS := OrderDetSymMatrixOfGroup(G);
  # Just drop the matrix.
  return ODS{[1..Length(ODS)]}{[1..3]};
end;

# We'll be referring to these lists a lot, so we may as well load them now.
if not IsBound(RodGroups) then
  # Avoid recalculating this every time I reload the file (most relevant for development).
  RodGroups := Immutable(List([1..75], RodGroupIT));
fi;
if not IsBound(RodGroupsOrderDetSym) then
  RodGroupsOrderDetSym := Immutable(List(RodGroups, OrderDetSymOfGroup));
fi;

RodGroupNumWithMatchingOps := function(G)
  local ODS;
  #% From a rod group, find the IT rod groups with the same representative operators up to chirality and axis.
  # Validate that we indeed have a Rod Group
  if not (IsAffineCrystGroupOnLeft(G) or IsAffineCrystGroupOnRight(G)) then
    Error("G is not an AffineCrystGroup");
  fi;
  if Length(TranslationBasis(G)) <> 1 or Length(TranslationBasis(G)[1]) <> 3 then
    Error("G is not a rod group");
  fi;
  # Check Order, Determinant, Symmomorphic
  ODS := OrderDetSymOfGroup(G);
  return Positions(RodGroupsOrderDetSym, ODS);
end;

# Step 2: Identify any glide planes or screw axes and match with the given rod
# groups. This would have to be done on a point group by point group basis.

# Alternative Step 2: find a transformation matrix that takes the (properly
# ordered) set of operators to the tabulated group. If that matrix exists,
# you have found the matching group.

# The non-distinct rod groups are in pairs.
# Either they have an axis with the translation axis and away from the translation axis, or they differ by chirality of a screw axis.

CrossProduct := function(a, b)
  # Assert that a and b are 3-vectors
  # Convert to row vectors if possible
  if not IsRowVector(a) then
    if IsMatrix(a) then
      if DimensionsMat(a) = [3,1] then
        a := a{[1..3]}[1];
      elif DimensionsMat(a) = [1,3] then
        a := a[1];
      fi;
    fi;
  fi;
  if not IsRowVector(b) then
    if IsMatrix(b) then
      if DimensionsMat(b) = [3,1] then
        b := b{[1..3]}[1];
      elif DimensionsMat(b) = [1,3] then
        b := b[1];
      fi;
    fi;
  fi;
  if not (IsRowVector(a) and IsRowVector(b)) then
    Error("CrossProduct only accepts vectors.");
  fi;
  if Length(a) <> 3 or Length(b) <> 3 then
    Error("CrossProduct requires 3-dimensional vectors.");
  fi;
  return [ a[2]*b[3] - a[3]*b[2],
           a[3]*b[1] - a[1]*b[3],
           a[1]*b[2] - a[2]*b[1] ];
end;

CrystGroupOnRightToStandardBasis := function(G)
  local T;
  T := Inverse(SymmetricInternalBasis(G));
  T := T * SignInt(Determinant(T)); # Preserve parity (assuming odd dimension)
  T := AugmentedMatrix(T, 0 * [1..Length(T)]);
  return G^T;
end;

FloorOfRational := r -> -CeilingOfRational(-r);

IsRodGroupEquivalentNC := function(G, Gref)
  local ODS, ODSref, screw, g, gref, grefs, op, v, x, wg, wgref;
  #% Compares two rod groups, determines if they're of the same type
  # Get to consistent action
  if IsAffineCrystGroupOnLeft(G) then
    G := TransposedMatrixGroup(G);
  fi;
  if IsAffineCrystGroupOnLeft(Gref) then
    Gref := TransposedMatrixGroup(Gref);
  fi;
  # Check trivial equality
  if G = Gref then
    return true;
  fi;
  # We'll assume that we're doing one of the pairs of groups which aren't unique in their ODS.
  # Orient the translation axes to be standard basis.
  G := CrystGroupOnRightToStandardBasis(G);
  Gref := CrystGroupOnRightToStandardBasis(Gref);
  # Get order, determinant, symmomorphic, operators
  ODS := OrderDetSymMatrixOfGroup(G);
  ODSref := OrderDetSymMatrixOfGroup(Gref);
  # Step 2: find operation to pair between G and Gref.
  # If there's a screw operation, pair highest order one with same translation
  # Otherwise, find the C2 operation.
  # Otherwise, find the mirror operation.
  if ForAny(ODS, x -> x[1] > 1 and x[2] = 1 and x[3] = false) then
    # Screw operation: det=1, nonsymmorphic, not pure translation.
    screw := true;
    # ODS is sorted in ascending Order (then Det, then Sym).
    # So we'll take the last screw operation
    g := Last(Filtered(ODS, x -> x[2] = 1 and x[3] = false))[4];
    # Now go over the possible two grefs
    grefs := Filtered(ODSref, x -> x[1] = Order(g{[1..3]}{[1..3]}) and x[2] = 1 and x[3] = false);
    grefs := List(grefs, x -> x[4]);
    # Get the intrinsic part of the screw operation
    wg := DecomposeMatrixTranslationOnRight(g)[1];
    for op in grefs do
      wgref := DecomposeMatrixTranslationOnRight(op)[1];
      # Check that they are the same modulo lattice vector
      # And since we're in standard basis, that's the first component
      if wg[1] - FloorOfRational(wg[1]) = wgref[1] - FloorOfRational(wgref[1]) then
        gref := op;
        break;
      fi;
    od;
    if not IsBound(gref) then
      Print("Could not find matching screw rotation.\n");
      Print("g: ");
      Display(g);
      Print("wg: ", DecomposeMatrixTranslationOnRight(g)[1], "\n");
      Print("grefs: ");
      Display(grefs);
      Print("wgrefs: ", List(grefs, g->DecomposeMatrixTranslationOnRight(g)[1]), "\n");
      return false;
    fi;
  elif ForAny(ODS, x -> x[1] = 2 and x[2] = 1) then
    # 2-fold rotation
    screw := false;
    # As C2 is its own inverse, there is only one C2.
    g := Filtered(ODS, x -> x[1] = 2 and x[2] = 1)[1][4];
    gref := Filtered(ODSref, x -> x[1] = 2 and x[2] = 1)[1][4];
  elif ForAny(ODS, x -> x[1] = 2 and x[2] = -1 and x[3] = true) then
    # Mirror operation
    screw := false;
    g := Filtered(ODS, x -> x[1] = 2 and x[2] = -1)[1][4];
    gref := Filtered(ODSref, x -> x[1] = 2 and x[2] = -1)[1][4];
  else
    Print("Could not find a screw operation, 2-fold rotation, or mirror operation in G.\n");
    return ODS{[1..Length(ODS)]}{[1..3]} = ODSref{[1..Length(ODSref)]}{[1..3]};
  fi;
  # Step 3a: If C2 or mirror, check if axis is parallel.
  # I.e. Does the operation leave translation vector unchanged?
  if not screw then
    v := ShallowCopy(TranslationBasis(G)[1]);
    Add(v, 0); # Convert 3D vector to augemented matrix form.
    # Check if rotation leaves the translation axis invariant
    # or rotation reflects the translation axis.
    # The groups are the same if they both do the same thing.
    return (v * Determinant(g) = v * g) = (v * Determinant(gref) = v * gref);
  else
    # Step 3b: If screw, check orientation/handedness.
    # To get the sense of the rotation, we need a vector orthogonal to the
    # translation basis.
    v := TranslationBasis(G)[1];
    if v <> [1,0,0] then
      x := [1,0,0];
    else
      x := [0,1,0];
    fi;
    x := x - x * v; # Subtract out component parallel to v
    # x is now orthogonal to v.
    return SignInt(CrossProduct(x, (x*g){[1..3]}) * v) = SignInt(CrossProduct(x, (x*gref){[1..3]}) * v);
  fi;
end;

# Finally, the function we've all been waiting for.
# The function which tells us which rod group G belongs to.
ITNumberOfRodGroup := function(G)
  local nums, num;
  nums := RodGroupNumWithMatchingOps(G);
  if Length(nums) = 1 then
    return nums[1];
  elif Length(nums) = 2 then
    for num in nums do
      if IsRodGroupEquivalentNC(G, RodGroupIT(num)) then
        return num;
      fi;
    od;
  fi;
  return fail;
end;

AxisOfOperator := function(M)
  local d, M2, det;
  # Gets the rotation axis or reflection normal of affine operator M.
  # Gives spurious results for plain inversion or the identity.
  d := Length(M) - 1;
  M2 := M{[1..d]}{[1..d]};
  det := Determinant(M2);
  # A rotation has determinant = 1.
  # Its axis is unchanged by rotation, so is an eigenvector of value 1.
  # A reflection has determinant = -1.
  # Its normal is reflected, so is an eigenvector of value -1.
  return GeneratorsOfVectorSpace(Eigenspaces(Rationals, M2)[Position(Eigenvalues(Rationals, M2), det)])[1];
end;

OrientationOfVector := function(v)
  # If v has one non-zero element, return index of that element. Else fail.
  if Length(v) - Length(Positions(v,0)) = 1 then
    return PositionProperty(v, x -> x <> 0);
  else
    return fail;
  fi;
end;

# We also want to identify the setting of the rod group.
SettingOfRodGroupNC := function(G, num)
  local T, c, gens, ax, a, M, b;
  # Group G, with IT number num.
  if (num >= 3 and num <= 22) then
    # Monoclinic and orthorhombic
    # Identify translation vector direction
    T := TranslationBasis(G);
    c := OrientationOfVector(T[1]);
    if c = fail then
      Print("TranslationBasis ",T," is not standard.\n");
      return fail;
    fi;
    # Case: the direction of the translation basis fully defines the setting.
    if num in [8,9,10,11,12,13,14,15,16,20,21] then
      if c = 3 then
        return "abc";
      elif c = 1 then
        return "cba";
      elif c = 2 then
        return "bca";
      else
        Error("You should not have gotten here. (return to continue)\n");
        return fail;
      fi;
    # Case 2: the representative symmetry ops only lie along the a axis
    elif num in [3..7] then
      # Get a generator which isn't translation, identity, or inversion
      gens := Filtered(GeneratorsOfGroup(G),
              x->x{[1..3]}{[1..3]} <> IdentityMat(3) and x{[1..3]}{[1..3]} <> - IdentityMat(3));
      ax := AxisOfOperator(gens[1]);
      a := OrientationOfVector(ax);
      if a = fail then
        Print("Axis of generator ",gens[1]," is ",ax,", which is non-standard.\n");
        return fail;
      fi;
    # Now 17,18,19,22
    elif num in [18,19] then
      # The 2-fold rotation (there's only one) points along vector a.
      # Order = 2, Determinant = 1, get the matrix.
      M := Filtered(OrderDetMatrixOfGroup(G), x -> x[1]=2 and x[2]=1)[1][3];
      ax := AxisOfOperator(M);
      a := OrientationOfVector(ax);
      if a = fail then
        Print("Axis of generator ",M," is ",ax,", which is non-standard.\n");
        return fail;
      fi;
    elif num in [17,22] then
      # The b-axis is the normal of the glide reflection.
      M := Filtered(OrderDetSymMatrixOfGroup(G), x -> x[1]=2 and x[2]=-1 and x[3]=false)[1][4];
      ax := AxisOfOperator(M);
      b := OrientationOfVector(ax);
      if b = fail then
        Print("Axis of generator ",M," is ",ax,", which is non-standard.\n");
        return fail;
      fi;
      if b = c then
        Error("Somehow the translation axis and other operator are aligned. Did you use the right number?\n");
        return fail;
      fi;
      a := Filtered([1..3], x -> not x in [b,c])[1];
    else
      ErrorNoReturn("You shouldn't be here. A number was forgotten.\n");
    fi;
    if a = c then
      Error("Somehow the translation axis and other operator are aligned. Did you use the right number?\n");
      return fail;
    fi;
    if c = 1 then
      if a = 2 then
        return "cab";
      else
        return "cba";
      fi;
    elif c = 2 then
      if a = 1 then
        return "acb";
      else
        return "bca";
      fi;
    else
      if a = 1 then
        return "abc";
      else
        return "bac";
      fi;
    fi;
  elif num in [35,37,38,41] then
    # Tetragonal
    Print("Tetragonal rod group settings not implemented.\n");
    return '1';
  elif num in [46..52] or num in [70,71,72,75] then
    # Trigonal, hexagonal
    Print("Trigonal and hexagonal rod group settings not implemented.\n");
    return '1';
  else
    # Cases with only one setting.
    return '1';
  fi;
  Error("Exhausted all options. You shouldn't be here. (return to continue)\n");
  return fail;
end;

# Now do the same thing for Layer Groups
if not IsBound(LayerGroups) then
  # Avoid recalculating this every time I reload the file (most relevant for development).
  LayerGroups := Immutable(List([1..80], LayerGroupIT));
fi;
if not IsBound(LayerGroupsOrderDetSym) then
  LayerGroupsOrderDetSym := Immutable(List(LayerGroups, OrderDetSymOfGroup));
fi;
if not IsBound(LayerGroupsWyckoffCount) then
  LayerGroupsWyckoffCount := Immutable(List(LayerGroups, x -> Length(WyckoffPositions(x))));
fi;

LayerGroupNumWithMatchingOps := function(G)
  local ODS;
  #% From a layer group, find the IT rod groups with the same representative operators up to chirality and axis.
  # Validate that we indeed have a Layer Group
  if not (IsAffineCrystGroupOnLeft(G) or IsAffineCrystGroupOnRight(G)) then
    Error("G is not an AffineCrystGroup");
  fi;
  if Length(TranslationBasis(G)) <> 2 or Length(TranslationBasis(G)[1]) <> 3 then
    Error("G is not a layer group");
  fi;
  # Check Order, Determinant, Symmomorphic
  ODS := OrderDetSymOfGroup(G);
  return Positions(LayerGroupsOrderDetSym, ODS);
end;

FilterLayerGroupNumWithMatchingWyckoff := function(G, nums)
  local nWyckoff;
  #% Group G and list of integers (IT numbers) nums. Returns subset of nums with same number of Wyckoff positions as G.
  nWyckoff := Length(WyckoffPositions(G));
  return Filtered(nums, i -> LayerGroupsWyckoffCount[i] = nWyckoff);
end;

TetragonalTrigonalLayerToStandardNC := function(G)
  local OD, W, T, t, T2;
  # Get order, determinant, operators (don't need expensive symmorphic info)
  OD := OrderDetMatrixOfGroup(G);
  # Check order information
  if ForAny(OD, x -> x[1] >= 3) then
    # Trigonal or square system. Find a 3-fold or 4-fold rotation
    W := OD[PositionProperty(OD, x -> (x[1] = 3) or (x[1] = 4))][3];
  else
    # If G isn't higher order,
    # don't need to do further basis transformation
    return G;
  fi;
  # Extract just the bare rotation
  W := W{[1..3]}{[1..3]} / Determinant(W);
  # Get the two translation vectors
  T := TranslationBasis(G);
  if (T[1] * W = T[2]) or (T[2] * W = T[1]) then
    # We're already in a suitable basis
    return G;
  fi;
  for t in T do
    T2 := [t, t * W, [0,0,1]];
    # Choose the basis vector which, when rotated, gives a determinant of 1.
    if Determinant(T2) = -1 then
      T2[2] := T2[2] * -1;
      break;
    fi;
    if Determinant(T2) = 1 then
      break;
    fi;
  od;
  if Determinant(T2) <> 1 then
    Print("Warning! Basis transformation has determinant ", Determinant(T2), "!\n");
  fi;
  # Apply the basis transformation
  T2 := AugmentedMatrix(T2, [0,0,0]);
  return G ^ Inverse(T2);
end;
IsGlideAlignedWithTwoFoldAxisNC := function(G)
  local ODS, W, axis, aligned, i, wg;
  #% Just for testing layer groups 28, 30, 32, or 34.
  # Get the matrices again
  ODS := OrderDetSymMatrixOfGroup(G);
  # Get the rotation axis of the 2-fold rotation.
  W := ODS[PositionProperty(ODS, x -> (x[1] = 2) and (x[2] = 1))][4];
  W := W{[1..3]}{[1..3]}; # Linear part
  # The axis is eigenvector with eigenvalue of 1
  axis := GeneratorsOfVectorSpace(Eigenspaces(Rationals, W)[Position(Eigenvalues(Rationals, W), 1)])[1];
  # Now we go over the non-symmorphic mirror ops, get their glide vectors,
  # and see if they align with axis.
  aligned := true;
  for i in PositionsProperty(ODS, x -> (x[1] = 2) and (x[2] = -1) and (x[3] = false)) do
    W := ODS[i][4]; # Mirror glide operator
    wg := DecomposeMatrixTranslationOnRight(W)[1]; # Glide vector
    # Test alignment with axis
    if (wg * axis) ^ 2 <> (wg * wg) * (axis * axis) then
      aligned := false;
      break;
    fi;
  od;
  return aligned;
end;

IsOrderDetSymMatrixInPlane := function(ODSM)
  local order, det, symmorphic, M, z;
  #% Checks whether an operation (given by list [order, determinant, symmorphic, matrix]) of a layer group in standard basis is in-plane and RightAction.
  # Unpack input
  order := ODSM[1];
  det := ODSM[2];
  symmorphic := ODSM[3];
  M := ODSM[4];
  # First trivial case: order is not 2.
  if order < 2 then
    # Identity or inversion
    return true;
  elif order > 2 then
    # Higher-order rotation. It's axis is necessarily out-of-plane.
    return false;
  fi;
  # Next case: screw rotation
  if det = 1 and not symmorphic then
    # As all glide vectors are necessarily in-plane, the screw roation axis is too
    return true;
  fi;
  # We should now be left with 2-fold rotations, mirror, and glide reflections.
  z := [0,0,1,0]; # Unit vector out-of-plane
  return (z + z*M)[3] = 0;
  # Applying the operation flips the z-component if true.
  # This means the mirror is in the plane (and its normal out of plane)
  # or the rotation axis is in-plane.
  # (I could define true to be if the normal is in-plane, but as long as I use
  # a consistent convention we're fine.)
end;

IsLayerGroupEquivalentNC := function(G, G2)
  local ODS, ODS2, func;
  #% Compares two layer groups, determines if they're of the same type, after filtering by operations and Wyckoff positions.
  # This algorithm works by categorising each operation as being in-plane or out-of-plane.
  # Get to consistent action
  if IsAffineCrystGroupOnLeft(G) then
    G := TransposedMatrixGroup(G);
  fi;
  if IsAffineCrystGroupOnLeft(G2) then
    G2 := TransposedMatrixGroup(G2);
  fi;
  # Check trivial equality
  if G = G2 then
    return true;
  fi;
  # We'll assume that we're doing one of the pairs of groups which aren't unique in their ODS.
  # Orient the translation axes to be standard basis.
  G := CrystGroupOnRightToStandardBasis(G);
  G2 := CrystGroupOnRightToStandardBasis(G2);
  # Except, in tetragonal (square) or trigonal/hexagonal systems, there is
  # an extra qualifier for "standard basis": that the basis rotates into itself.
  #G := TetragonalTrigonalLayerToStandardNC(G);
  #G2 := TetragonalTrigonalLayerToStandardNC(G2);
  # Get order, determinant, symmomorphic, operators
  ODS := OrderDetSymMatrixOfGroup(G);
  ODS2 := OrderDetSymMatrixOfGroup(G2);
  # Check orientations of operations
  # Drop matrices
  func := function(lst)
    Add(lst, IsOrderDetSymMatrixInPlane(lst));
    Remove(lst, 4);
  end;
  Perform(ODS, func);
  Perform(ODS2, func);
  # Sort
  Sort(ODS);
  Sort(ODS2);
  # Compare
  #Print("ODS: ");
  #Display(ODS);
  #Print("ODS2: ");
  #Display(ODS2);
  return ODS = ODS2;
end;

LayerGroupNum31Or36NC := function(G)
  local ODSM, M, x, v1, v2, T;
  # Standardise G
  if IsAffineCrystGroupOnLeft(G) then
    G := TransposedMatrixGroup(G);
  fi;
  G := CrystGroupOnRightToStandardBasis(G);
  # Get the operations and their details
  ODSM := OrderDetSymMatrixOfGroup(G);
  # Step one: find lattice vector along the glide vector
  # Find the glide operation
  for x in ODSM do
    if x[1] = 2 and x[2] = -1 and not x[3] then
      M := x[4];
      break;
    fi;
  od;
  if not IsBound(M) then
    Print("Could not find glide operation in group.\n");
    return fail;
  fi;
  # Now get the glide vector (or, specifically, double the glide vector), which shall be a new basis vector
  v1 := (M * M)[4]{[1..3]};
  # Step two: find lattice vector along the 2-fold rotation axis
  Unbind(M);
  for x in ODSM do
    if x[1] = 2 and x[2] = 1 then
      M := x[4];
      break;
    fi;
  od;
  if not IsBound(M) then
    Print("Could not find rotation operation in group.\n");
    return fail;
  fi;
  # Find the axis of the 2-fold rotation.
  v2 := AxisOfOperator(M);
  if CrossProduct(v1, v2) = [0,0,0] then
    # If v1 and v2 are parallel, then
    # This can only happen in group cm2e.
    return 36;
  fi;
  # Now we normalise v2 to be the smallest lattice vector.
  # Because we are in standard basis, and everything so far should be integers, we just divide by
  # the greatest common divisor
  if not ForAll(v2, IsInt) then
    # Ensure v2 is integers by multiplying out denominators.
    v2 := v2 * MaximumList(List(v2, DenominatorRat));
  fi;
  v2 := v2 / Gcd(Integers, v2);
  # Now we change basis
  M := IdentityMat(4);
  M[1]{[1..3]} := v1;
  M[2]{[1..3]} := v2;
  M := Inverse(M);
  G := G ^ M;
  # Now we check the *primitive* translation basis
  T := TranslationBasis(G);
  # If this basis has fractional components, then we have a centering lattice
  if 1/2 in Flat(T) then
    return 36;
  else
    return 31;
  fi;
end;

ITNumberOfLayerGroup := function(G)
  local nums, num;
  # First, get the layer groups with have matching operation types.
  nums := LayerGroupNumWithMatchingOps(G);
  # If we have multiple operations, also filter by number of Wyckoff positions.
  if Length(nums) > 1 then
    nums := FilterLayerGroupNumWithMatchingWyckoff(G, nums);
  fi;
  if Length(nums) = 1 then
    # If we only have one group left, it must be that one
    return nums[1];
  elif Length(nums) > 1 then
    # Otherwise, we need to narrow it down further
    if nums = [31,36] then
      # This one, the comparison between 31 and 36, is a special case
      # Handle it separately
      return LayerGroupNum31Or36NC(G);
    fi;
    # Otherwise, the remaining sets of groups can be distinguished by checking whether
    # their respective operations are in-plane or out-of-plane.
    for num in nums do
      # Could use LayerGroupIT(num) instead of LayerGroups[num], but
      # I preloaded this for a reason (specifically performance).
      if IsLayerGroupEquivalentNC(G, LayerGroups[num]) then
        return num;
      fi;
    od;
  fi;
  return fail;
end;

# We also want to identify the setting of the layer group
SettingOfLayerGroupNC := function(G, num)
  local ops, M, T, OD, t, a, ax, b, Tdir, wg;
  # Layer group G with IT number num
  if IsAffineCrystGroupOnLeft(G) then
    # Set to right-action
    G := TransposedMatrixGroup(G);
  fi;
  if num = 5 or num = 7 then
    # p11a or p112/a
    # Get the mirror glide vector
    OD := OrderDetMatrixOfGroup(G);
    M := Filtered(OD, x -> x[1] = 2 and x[2] = -1)[1][3];
    wg := DecomposeMatrixTranslationOnRight(M)[1];
    # Convert wg to standard basis
    wg := wg * Inverse(SymmetricInternalBasis(G));
    # Normalise wg to unit cell
    wg := List(wg, x -> x - FloorOfRational(x));
    # Compare with translation basis
    if wg = [1/2,0,0] then
      return '1';
    elif wg = [1/2,1/2,0] then
      return '2';
    elif wg = [0,1/2,0] then
      return '3';
    else
      Print("Group L",num," has non-standard glide vector ",wg," in basis ",TranslationBasis(G),"\n");
      return fail;
    fi;
  elif num >= 8 and num <= 48 then
    # Rectangular systems
    if num in [19,21,22,23,25,26,37,39,44,46,47,48] then
      # These rectangular systems have only one setting.
      return '1';
    elif num in [10,13,18,35,36] then
      # Groups with centering, so require special care.
      # International Tables choose a standard basis which is non-primitive,
      # but GAP typically translates into a primitive standard basis.
      # I want to be able to catch both cases.
      T := TranslationBasis(G);
      Tdir := List(T, OrientationOfVector);
      # Get the operator of interest
      if num in [10,13] then
        # c211 and cm11. Only one non-translation operator.
        M := Filtered(GeneratorsOfGroup(G),
             x -> x{[1..3]}{[1..3]} <> IdentityMat(3))[1];
      else
        # c2/m11, cm2m, and cm2e. We want the 2-fold rotation.
        OD := OrderDetMatrixOfGroup(G);
        M := Filtered(OD, x -> x[1] = 2 and x[2] = 1)[1][3];
      fi;
      ax := AxisOfOperator(M);
      # Now we determine the direction of this axis
      t := OrientationOfVector(ax);
      if t = fail then
        # This may occur if our basis is primitive.
        if ForAll(Tdir, IsInt) then
          # The basis is primitive, so correct for it.
          # The lattice vector which is orthogonal to ax will be another eigenvector of M,
          # with eigenvalue of -Determinant(M).
          t := PositionProperty(T, v -> v * M{[1..3]}{[1..3]} = -Determinant(M)*v);
        fi;
        if num in [10,13,18] then
          b := t;
        else
          a := t;
        fi;
      elif num in [10,13,18] then
        a := t;
      else
        b := t;
      fi;
    elif num <= 18 then
      # Only one axis has symmetry operations
      ops := Filtered(GeneratorsOfGroup(G),
             x -> x{[1..3]}{[1..3]} <> IdentityMat(3) and x{[1..3]}{[1..3]} <> - IdentityMat(3));
      ax := AxisOfOperator(ops[1]);
      a := OrientationOfVector(ax);
    elif num = 20 then
      # p2_122
      # Get the screw rotation.
      OD := OrderDetSymMatrixOfGroup(G);
      M := Filtered(OD, x -> not x[3])[1][4];
      # Get its screw vector
      ax := DecomposeMatrixTranslationOnRight(M)[1];
      # Normalise wg to unit cell
      ax := List(ax, x -> x - FloorOfRational(x));
      a := OrientationOfVector(ax);
    elif num in [24, 31, 38, 40,41] or num in [28,29,30] then
      # Groups with glide(s) that are just in one direction.
      # Get the glide operation
      OD := OrderDetSymMatrixOfGroup(G);
      M := Filtered(OD, x -> (not x[3]) and x[2] = -1)[1][4];
      # Get its glide vector
      ax := DecomposeMatrixTranslationOnRight(M)[1];
      # Normalise wg to unit cell
      ax := List(ax, x -> x - FloorOfRational(x));
      if num in [28,29,30] then
        # The glide is in the b direction
        b := OrientationOfVector(ax);
      else
        # The glide is in the a direction
        a := OrientationOfVector(ax);
        # This might fail for L36, cm2a. In fact, it likely will.
      fi;
    elif num in [27,32,33,34] then
      # These groups have just one 2-fold rotation, along the b direction
      # Find the 2-fold rotation
      OD := OrderDetMatrixOfGroup(G);
      M := Filtered(OD, x -> x[1] = 2 and x[2] = 1)[1][3];
      ax := AxisOfOperator(M);
      b := OrientationOfVector(ax);
    elif num in [42,43] then
      # 2-fold screw is in b direction
      OD := OrderDetSymMatrixOfGroup(G);
      M := Filtered(OD, x -> x[1] = 2 and x[2] = 1 and x[3] = false)[1][4];
      ax := AxisOfOperator(M);
      b := OrientationOfVector(ax);
    elif num = 45 then
      # pbma. m points in b direction
      OD := OrderDetSymMatrixOfGroup(G);
      M := Filtered(OD, x -> x[1] = 2 and x[2] = -1 and x[3] = true)[1][4];
      ax := AxisOfOperator(M);
      b := OrientationOfVector(ax);
    else
      Error("Did not catch case for number ",num,". You should not be here. (return to continue)\n");
      return fail;
    fi;
    if IsBound(b) then
      if b = 1 then
        a := 2;
      elif b = 2 then
        a := 1;
      else
        a := fail;
      fi;
    fi;
    if a = 1 then
      return 'a';
    elif a = 2 then
      return 'b';
    else
      Print("Group L",num," has non-standard axis ",ax,"\n");
      return fail;
    fi;
  elif num in [52, 62, 64] then
    # Two origin choices
    # Check we are in standard basis
    T := TranslationBasis(G);
    if not ([1,0,0] in T and [0,1,0] in T) then
      Print("Layer group has non-standard translation basis ",T,"\n");
      return fail;
    fi;
    # Grab the inversion operator
    OD := OrderDetMatrixOfGroup(G);
    M := Filtered(OD, x -> x[1] = -1)[1][3];
    t := M[4]{[1..3]};
    # Normalise to unit cell
    t := List(t, x -> x - FloorOfRational(x));
    if t = [1/2, 1/2, 0] then
      return '1';
    elif t = [0, 0, 0] then
      return '2';
    else
      Print("Group L",num," has non-standard inversion translation ",t,"\n");
      return fail;
    fi;
  else
    # Cases with only one setting
    return '1';
  fi;
  Error("Exhausted all options. You shouldn't be here. (return to continue)\n");
  return fail;
end;
