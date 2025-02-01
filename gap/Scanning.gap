#!/usr/bin/env gap
# Copyright 2024 Bernard Field

LoadPackage("NumericalSgps",false); # Includes CeilingOfRational method.
#Read("RepresentativeSymmetryOps.gap");
Read("Identify_Group.gap"); # RepresentativeSymmetryOp.gap is read in this one

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

VectorNormSquared := function(v)
  if not IsRowVector(v) then
    Error("v is not a vector.");
  fi;
  return Sum(List(v, x->x^2));
end;

VectorLengthRatio := function(v1, v2)
  local Rsq, num, den;
  #% ratio of |v2|/|v1|. Problems if irrational.
  # Square of the ratio
  Rsq := VectorNormSquared(v2) / VectorNormSquared(v1);
  if IsInt(Rsq) then
    if IsSquareInt(Rsq) then
      return RootInt(Rsq);
    else
      Print("Ratio squared, ", Rsq, ", is not a square number! Returning a float.\n");
      return Rsq^0.5;
    fi;
  else
    # Rsq is rational, so do the numerator and denominator separately.
    num := NumeratorRat(Rsq);
    den := DenominatorRat(Rsq);
    if IsSquareInt(num) and IsSquareInt(den) then
      return RootInt(num) / RootInt(den);
    else
      Print("Ratio squared, ", Rsq, ", is not a square number! Returning a float.\n");
      return Rsq^0.5;
    fi;
  fi;
end;

StandardOriginTransformationOfGroupByInversion := function(G)
  local d, g, T, t, s;
  # Per International Tables, if it's centrosymmetric, put origin on an
  # inversion center
  d := DimensionOfMatrixGroup(G) - 1; # Dimension
  # Find the inversion operator
  for g in RepresentativeSymmetryOps(G) do;
    if g{[1..d]}{[1..d]} = -1*IdentityMat(d) then;
      # Shift the origin to cancel out the inversion's shift
      if IsAffineCrystGroupOnLeft(G) then
        t := g{[1..d]}[d+1];
      else
        t := g[d+1]{[1..d]};
      fi;
      # If t is a lattice vector, we don't want that.
      T := TranslationBasis(G);
      if IntSolutionMat(T, t) <> fail then
        return IdentityMat(d+1);
      else
        # Shift t to be in the home unit cell.
        s := SolutionMat(T,t);
        # There's a chance that the origin is currently outside the subperiodic
        # group's basis. In that case, not worth doing anything fancy.
        # But for nicer cases, shifting to home unit cell allows for better
        # comparisons between differently-expressed-but-otherwise-equal groups.
        if s <> fail then
          s := List(SolutionMat(T,t), FloorOfRational);
          t := t - s * T;
        fi;
        return TranslationAffineMatrix(-t/2);
      fi;
    fi;
  od;
  Error("Did not find inversion center in when searching of inversion center for standard origin.");
  return IdentityMat(d+1);
end;

StandardOriginTransformationOfGroupByWyckoff := function(G)
  local d, W, wy, maxwy, wysize, maxsize, needtranslate, needtranslatetmp, t, M, x;
  # For non-centrosymmetric groups, the origin is at a point of a highest symmetry site.
  d := DimensionOfMatrixGroup(G) - 1;
  # First, get our symmetry sites, the Wyckoff Positions
  W := WyckoffPositions(G);
  # Now find those with the greatest symmetry (i.e. those with the biggest
  # stabiliser/point group).
  maxsize := 0;
  maxwy := W[1];
  needtranslate := true;
  for wy in W do;
    wysize := Size(WyckoffStabilizer(wy));
    if (wysize > maxsize) or (wysize = maxsize and needtranslate) then
      # If we have a new biggest, accept it unconditionally.
      # Otherwise, if we have a tie, only accept it if our current best
      # needs a translation (in case we can find one which doesn't).
      
      # Check whether a translation will be necessary
      # This is done by checking whether this Wyckoff position is equivalent
      # to a similar Wyckoff position which intersects the origin
      needtranslatetmp := WyckoffPositionObject(rec(basis:=WyckoffBasis(wy), class:=1, spaceGroup := WyckoffSpaceGroup(wy), translation := [1..d]*0)) <> wy;
      # If this new Wyckoff position is of the same size and also needs
      # translating, then only accept it if the vector is shorter.
      # This shortest-if-all-else-equal condition ensures we get the same answer
      # even when the positions are arranged differently.
      if (wysize > maxsize) or not needtranslatetmp or WyckoffTranslation(wy)*WyckoffTranslation(wy) < WyckoffTranslation(maxwy)*WyckoffTranslation(maxwy) then
        maxwy := wy;
        maxsize := wysize;
        needtranslate := needtranslatetmp;
      fi;
    fi;
  od;
  if verbose then; Print("Highest symmetry site: ", maxwy, "\n"); fi;
  # Now, translate such that this Wyckoff position intersects the new origin.
  if needtranslate then
    # We need to be careful. It isn't always as simple as subtracting off
    # translation, because it might have an origin not on the rod.
    # This is not a concern for space groups, or if WyckoffBasis is a point,
    # or if WyckoffTranslation lies in the span of TranslationBasis(G).
    # The latter is tested with the Rouche-Capelli theorem.
    if IsSpaceGroup(G) or Length(WyckoffBasis(maxwy)) = 0 or
      (RankMatrix(TransposedMat(TranslationBasis(G))) = RankMatrix(TransposedMat(Concatenation(TranslationBasis(G), [WyckoffTranslation(maxwy)]))))
      then
      t := WyckoffTranslation(maxwy);
    elif Length(TranslationBasis(G)) = 1 or Length(WyckoffBasis(maxwy)) = 1 then
      # G is a rod group or Wyckoff position is a line. We have a unique solution.
      # Solve wantedTranslation = currentTranslation + WyckoffBasis*x2 = TranslationBasis * x1
      # (x1|x2) = (TranslationBasis|-WyckoffBasis)^-1 * currentTranslation
      M := TransposedMat(Concatenation(TranslationBasis(G),-WyckoffBasis(maxwy)));
      # Since M may be rectangular, we instead compute the pseudoinverse
      # to get the left inverse.
      x := (TransposedMat(M)*M)^-1*TransposedMat(M) * WyckoffTranslation(maxwy);
      t := TransposedMat(TranslationBasis(G)) * x{[1..Length(TranslationBasis(G))]};
    else
      # We have a non-unique solution. This could be problematic.
      # The proper thing to do is to find the intersection of several
      # Wyckoff positions (although isn't that just a higher symmetry
      # Wyckoff position?). I'm not sure what a good procedure for that is.
      # However, I don't anticipate creating off-axis groups.
      t := WyckoffTranslation(maxwy);
      Error("Situation I can't handle in StandardOriginTransformationOfGroup, a subperiodic group where the WyckoffPosition has an off-axis WyckoffTranslation and no unique choice of way to get it back on-axis.\n");
    fi;
    return TranslationAffineMatrix(-t);
  fi;
  # needtranslate = false, no translation necessary
  return IdentityMat(d+1);
end;

StandardOriginTransformationOfLayerGroup48 := function(G)
  local OD, ops, axs, adirs, P, action, x1, x2, T, P2, idx, tdirs, t, i, a;
  # cmme is a special case: it has two inversion centers with distinct settings.
  # The "standard" origin lies on mx.
  # Of course, this only matters if we're in something resembling a standard setting.
  # Grab the two mirror operators.
  OD := OrderDetSymMatrixOfGroup(G);
  ops := Filtered(OD, x -> x[1] = 2 and x[2] = -1 and x[3]);
  ops := ops{[1..Length(ops)]}[4];
  if Length(ops) <> 2 then
    Error("We don't have 2 mirror operations. I don't think this group is L48.");
  fi;
  # Get the directions of the mirror normals.
  axs := List(ops, AxisOfOperator);
  adirs := List(axs, OrientationOfVector);
  # Go through the cases
  if ForAll(adirs, x -> x = fail) or ForAny(adirs, x -> x = 3) then
    # If we are in nothing resembling a standard setting, do the regular thing.
    return StandardOriginTransformationOfGroupByInversion(G);
  else
    # We are either in a conventional or standard-primitive basis.
    if IsAffineCrystGroupOnLeft(G) then
      action := LeftAction;
    else
      action := RightAction;
    fi;
    # Find where an inversion center is.
    P := StandardOriginTransformationOfGroupByInversion(G);
    # Remember, this transformation is what's needed to move the origin to 0.
    # So we want the negative of it.
    if action = LeftAction then
      x1 := -P{[1..3]}[4];
    else
      x1 := -P[4]{[1..3]};
    fi;
    # Determine which mirror operation x1 lies on (modulo lattice vector).
    T := TranslationBasis(G);
    for i in [1..Length(ops)] do
      # Have the mirror operation act on this point.
      if action = LeftAction then
        P2 := ops[i] * Inverse(P);
        x2 := P2{[1..3]}[4];
      else
        P2 := Inverse(P) * ops[i];
        x2 := P2[4]{[1..3]};
      fi;
      # If this point is left invariant modulo lattice vector, then we have our operation.
      if IntSolutionMat(T, x1 - x2) <> fail then
        idx := i;
        break;
      fi;
    od;
    if not IsBound(idx) then
      Error("Somehow L48 doesn't have an inversion center on a mirror plane...");
      return P;
    fi;
    # Now we determine if we want this inversion center or the other one.
    if adirs[idx] = 1 or (adirs[idx] = fail and adirs[(idx mod 2)+1] = 2) then
      # the operation of interest is mirror plane along the x axis.
      # This is either shown by this operation being along the first axis,
      # or the OTHER operation being along the y axis.
      return P;
    else
      # We want the OTHER inversion center.
      # Shift by half the centering vector.
      if ForAll(adirs, IsInt) then
        # We're in a conventional basis. The centering vector is the member of
        # T which isn't properly oriented.
        tdirs := List(T, OrientationOfVector);
        t := T[Position(tdirs, fail)];
      else
        # We're in a primitive basis.
        # Find the lattice vector NOT aligned with the aligned mirror axis.
        a := adirs[PositionProperty(adirs, IsInt)];
        t := T[PositionProperty(T, x -> OrientationOfVector(x) <> a)];
      fi;
      x1 := x1 + t/2;
      P := IdentityMat(4);
      if action = LeftAction then
        P{[1..3]}[4] := -x1;
      else
        P[4]{[1..3]} := -x1;
      fi;
      return P;
    fi;
  fi;
end;

StandardOriginTransformationOfRodGroup14 := function(G)
  local P, basis, c, ops, op, 2x, 2xS, Q, t, t2;
  # The origin sits on the 2-fold axis pointing along 'x'.
  # This is, of course, only meaningful if G is in a Cartesian basis.
  # Get one of the candidate origins.
  P := StandardOriginTransformationOfGroupByWyckoff(G);
  # Check the basis.
  c := TranslationBasis(G)[1];
  if OrientationOfVector(c) <> 3 then
    # The c-vector should be this one
    # In principle, other settings exist.
    # If I had to, I could transform G into a setting with c pointing in this direction
    # or do a setting-dependent check.
    # But I don't think I need that for my purposes.
    return P;
  fi;
  basis := SymmetricInternalBasis(G);
  if 1 in List(basis, OrientationOfVector) then
    # We have a well defined x axis.
    # Find the 2x operator.
    ops := RepresentativeSymmetryOps(G);
    for op in ops do
      if Order(op) > 1
        and ((IsAffineCrystGroupOnRight(G) and [1,0,0,0] * op = [1,0,0,0]) 
        or (IsAffineCrystGroupOnLeft(G) and op * [1,0,0,0] = [1,0,0,0])) then
        2x := op;
        break;
      fi;
    od;
    if not IsBound(2x) then
      Error("Basis of R14 includes [1,0,0] but unable to find 2x operator! (Return to continue with previously-found origin.");
      return P;
    fi;
    # We now check if the origin in P is on 2x or not.
    if IsAffineCrystGroupOnRight(G) then
      2xS := Inverse(P) * 2x * P;
      t := 2xS[4]{[1..3]};
    else
      2xS := P * 2x * Inverse(P);
      t := 2xS{[1..3]}[4];
    fi;
    if IsInt(t[3]/c[3]) then
      # We have the correct origin.
      return P;
    else
      # We have the wrong origin, by c/4.
      if IsAffineCrystGroupOnRight(G) then
        t := P[4]{[1..3]};
      else
        t := P{[1..3]}[4];
      fi;
      # It's +/- c/4, so choose one which takes us closer to zero.
      if (t[3] < 0) = (c[3] > 0) then
        t2 := c[3]/4;
      else
        t2 := -c[3]/4;
      fi;
      # Create a new transformation matrix.
      Q := List(P, ShallowCopy); # Clone P.
      if IsAffineCrystGroupOnRight(G) then
        Q[4][3] := Q[4][3] + t2;
      else
        Q[3][4] := Q[3][4] + t2;
      fi;
      return Q;
    fi;
  else
    return P;
  fi;
end;

StandardOriginTransformationOfRodGroup19 := function(G)
  local P, 2x, 2xS, Q, c, t, t2, OD, t3, i;
  # The standard origin of R19 p2cm sits on the two-fold rotation axis, not the mirror
  # Get the default answer. Then we'll shift by c/4 if it's off.
  P := StandardOriginTransformationOfGroupByWyckoff(G);
  # Find the two-fold rotation
  OD := OrderDetMatrixOfGroup(G);
  for i in OD do
    if i[1] = 2 and i[2] = 1 then
      2x := i[3];
      break;
    fi;
  od;
  if not IsBound(2x) then
    Error("Group believed to be R19 did not contain a 2-fold rotation! (Return to use old origin.)");
    return P;
  fi;
  # Get the translation component
  if IsAffineCrystGroupOnRight(G) then
    2xS := Inverse(P) * 2x * P;
    t := 2x[4]{[1..3]};
  else
    2xS := P * 2x * Inverse(P);
    t := 2x{[1..3]}[4];
  fi;
  c := TranslationBasis(G)[1];
  if IsInt(VectorLengthRatio(c, t)) then
    # t is an integer multiple of c
    # We have the right origin.
    return P;
  else
    # We're off by c/4.
    if IsAffineCrystGroupOnRight(G) then
      t := P[4]{[1..3]};
    else
      t := P{[1..3]}[4];
    fi;
    # It's +/- c/4, so choose one which takes us closer to zero.
    t2 := t + c/4;
    t3 := t - c/4;
    if Sum(t2, x -> x^2) >= Sum(t3, x -> x^2) then
      # Choose the smaller vector
      # Or the one that's more likely to be negative if equal
      # (because origin is negative of this vector).
      t2 := t3;
    fi;
    # Copy the matrix.
    Q := List(P, ShallowCopy);
    if IsAffineCrystGroupOnRight(G) then
      Q[4]{[1..3]} := t2;
    else
      Q{[1..3]}[4] := t2;
    fi;
    return Q;
  fi;
end;

StandardOriginTransformationOfGroup := function(G)
  local d, nums;
  #% Gives the transformation matrix which shifts a group to have a standard origin (as defined by the International Tables).
  d := DimensionOfMatrixGroup(G) - 1;
  # Special case: layer group 47 cmmm has two distinct inversion
  # centers. The standard choice is on a highest symmetry Wyckoff position.
  if d = 3 and Length(TranslationBasis(G)) = 2 then
    # Layer group
    nums := LayerGroupNumWithMatchingOps(G);
    if 47 in nums then
      nums := FilterLayerGroupNumWithMatchingWyckoff(G, nums);
      if nums = [47] then
        return StandardOriginTransformationOfGroupByWyckoff(G);
      fi;
    elif nums = [48] then
      # Another special case, L48 cmme has two inversion centers with distinct settings.
      return StandardOriginTransformationOfLayerGroup48(G);
    fi;
  # Special case: rod groups 14 p222_1 and 19 p2cm have two possible origin choices
  elif d = 3 and Length(TranslationBasis(G)) = 1 then
    # Rod group
    nums := RodGroupNumWithMatchingOps(G);
    if 14 in nums then
      return StandardOriginTransformationOfRodGroup14(G);
    elif 19 in nums then
      return StandardOriginTransformationOfRodGroup19(G);
    fi;
  fi;
  if -1*IdentityMat(d) in PointGroup(G) then
    # Per International Tables, if it's centrosymmetric, put origin on an
    # inversion center
    return StandardOriginTransformationOfGroupByInversion(G);
  else
    # Otherwise, find the origin using the Wyckoff positions.
    return StandardOriginTransformationOfGroupByWyckoff(G);
  fi;
end;

OriginLocationOfGroup := function(G)
  local M, v, T;
  # Get where the origin is
  M := StandardOriginTransformationOfGroup(G);
  if IsAffineCrystGroupOnLeft(G) then
    v := -M{[1..3]}[4];
  else
    v := -M[4]{[1..3]};
  fi;
  if SolutionMat(TranslationBasis(G), v) <> fail then
    # Wrap to positive fractional coordinates.
    # But don't wrap out-of-plane components.
    T := SymmetricInternalBasis(G);
    v := v * Inverse(T);
    v := List(v, x -> x - FloorOfRational(x));
    v := v * T;
  fi;
  return v;
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
