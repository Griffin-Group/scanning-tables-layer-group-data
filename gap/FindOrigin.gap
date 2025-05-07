#!/usr/bin/env gap
# Copyright 2025 Bernard Field

Read("Identify_Group.gap");

verbose := false;

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
  local d, W, wy, maxwy, wysize, maxsize, needtranslate, needtranslatetmp, t, M, x, accept;
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
    if (wysize >= maxsize) then
      # We have a new contender. If it's bigger, we'll accept it unconditionally.
      # But otherwise, we have a few criteria.
      # There are some Wyckoff positions that are preferred over others.
      # E.g. 4 > 222, or 2 > m, even if they both have the same size.
      # Otherwise, if the current position requires a translation but this new
      # position can do a smaller translation, do it.

      # Check whether a translation will be necessary
      # This is done by checking whether this Wyckoff position is equivalent
      # to a similar Wyckoff position which intersects the origin
      needtranslatetmp := WyckoffPositionObject(rec(basis:=WyckoffBasis(wy), class:=1, spaceGroup := WyckoffSpaceGroup(wy), translation := [1..d]*0)) <> wy;

      accept := false;
      if wysize > maxsize then
        # The new position has a greater size.
        accept := true;
      elif Maximum(List(WyckoffStabilizer(wy),Order)) > Maximum(List(WyckoffStabilizer(maxwy),Order)) then
        # The new position is on a higher-order rotation axis.
        accept := true;
      elif Maximum(List(WyckoffStabilizer(wy),Order)) < Maximum(List(WyckoffStabilizer(maxwy),Order)) then
        accept := false;
      elif Length(WyckoffBasis(wy)) < Length(WyckoffBasis(maxwy)) then
        # The position is on a rotation axis rather than a mirror plane.
        accept := true;
      elif Length(WyckoffBasis(wy)) > Length(WyckoffBasis(maxwy)) then
        accept := false;
      elif needtranslate then
        # If this new Wyckoff position is of the same size and also needs
        # translating, then only accept it if the vector is shorter.
        # This shortest-if-all-else-equal condition ensures we get the same answer
        # even when the positions are arranged differently.
        if not needtranslatetmp or WyckoffTranslation(wy)*WyckoffTranslation(wy) < WyckoffTranslation(maxwy)*WyckoffTranslation(maxwy) then
          accept := true;
        fi;
      fi;
      
      if accept then
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
    t := 2xS[4]{[1..3]};
  else
    2xS := P * 2x * Inverse(P);
    t := 2xS{[1..3]}[4];
  fi;
  c := TranslationBasis(G)[1];
  if IsInt(VectorLengthRatio(c, t)) then
    # t is an integer multiple of c
    # We have the right origin.
    # Or, at least, a right origin. But the origin is defined modulo c/2.
    # So we can potentially simplify.
    if IsAffineCrystGroupOnRight(G) then
      t := P[4]{[1..3]};
    else
      t := P{[1..3]}[4];
    fi;
    t2 := t - c/2;
    if Sum(t, x->x^2) > Sum(t2, x->x^2) then
      return TranslationAffineMatrix(t2);
    fi;
    t2 := t + c/2;
    if Sum(t, x->x^2) > Sum(t2, x->x^2) then
      return TranslationAffineMatrix(t2);
    fi;
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

StandardOriginTransformationOfLayerGroup4FoldTwoChoice := function(G)
  # L52, L62, and L64, with our origin choice 2 has the origin on an inversion center.
  # But the two inversion centres are subtly different.
  # The one we want is 1/4,1/4,0 from the Wyckoff position with 4, not -4.
  local W, wy, mywy, t, T, adirs, t1, Q, candidates, idx;
  
  # The procedure only makes sense if we are in something resembling an x-y translation basis.
  # If it is rotated, then we cannot distinguish the two settings/origins.
  T := TranslationBasis(G);
  adirs := List(T, OrientationOfVector);
  if ForAny(adirs, x -> x = fail) or ForAny(adirs, x -> x = 3) then
    # If we are in nothing resembling a standard setting, do the regular thing.
    return StandardOriginTransformationOfGroupByInversion(G);
  fi;

  # Find the Wyckoff position with 4.
  W := WyckoffPositions(G);
  for wy in W do
    if [4,1] in List(WyckoffStabilizer(wy), OrderDetOfMatrix) then
      mywy := wy;
      break;
    fi;
  od;
  if not IsBound(mywy) then
    # If G is indeed L62 then this should not happen.
    Error("Could not find Wyckoff position with 4 in claimed layer group 52/62/64! (return to use default inversion centre)");
    return StandardOriginTransformationOfGroupByInversion(G);
  fi;
  # Shift by 1/4,1/4 (or equivalent);
  t := -WyckoffTranslation(mywy);
  candidates := [t-T[1]/4-T[2]/4, t+T[1]/4+T[2]/4, t-3*T[1]/4-3*T[2]/4, t+3*T[1]/4+3*T[2]/4];
  # Find the shortest vector satisfying this.
  idx := PositionMinimum(candidates, x -> Sum(x, y->y^2));
  t1 := candidates[idx];
  Q := TranslationAffineMatrix(t1);
  return Q;
end;

StandardOriginTransformationOfGroupOn2z := function(G)
  local T, ST, z, W, needtranslate, needtranslatetmp, wy, mywy, t, M, x, d;
  # The standard origin of L20 p2_122 sits on the 2z axis, not the 2y axis.
  # L24 also has the origin sitting on 2z axes.
  # So we'll do the Wyckoff position thing, but instead on just the z-direction.

  # Find the z direction.
  # Fortunately for us, the 2z axis is well-defined as the out-of-plane direction.
  # Check the translation basis
  T := TranslationBasis(G);
  ST := SymmetricInternalBasis(G);
  z := Filtered(ST, x -> not x in T)[1];

  d := DimensionOfMatrixGroup(G) - 1;
  # Find the Wyckoff positions in a basis in the z direction.
  W := WyckoffPositions(G);
  needtranslate := true;
  for wy in W do
    # This will be our 2z.
    if WyckoffBasis(wy) = [z] then
      # Check whether a translation will be necessary
      # This is done by checking whether this Wyckoff position is equivalent
      # to a similar Wyckoff position which intersects the origin
      needtranslatetmp := WyckoffPositionObject(rec(basis:=WyckoffBasis(wy), class:=1, spaceGroup := WyckoffSpaceGroup(wy), translation := [1..d]*0)) <> wy;
      if not IsBound(mywy) then
        mywy := wy;
        needtranslate := needtranslatetmp;
      elif needtranslate then
        # If this new Wyckoff position is of the same size and also needs
        # translating, then only accept it if the vector is shorter.
        # This shortest-if-all-else-equal condition ensures we get the same answer
        # even when the positions are arranged differently.
        if not needtranslatetmp or WyckoffTranslation(wy)*WyckoffTranslation(wy) < WyckoffTranslation(mywy)*WyckoffTranslation(mywy) then
          mywy := wy;
          needtranslate := needtranslatetmp;
        fi;
      fi;
    fi;
  od;
  if not IsBound(mywy) then
    ErrorNoReturn("Could not find Wyckoff position with basis ", [z], " in group with Wyckoffs=", W);
  fi;
  # Now, translate such that this Wyckoff position intersects the new origin.
  if needtranslate then
    # We need to be careful. It isn't always as simple as subtracting off
    # translation, because it might have an origin not on the rod.
    # This is not a concern for space groups, or if WyckoffBasis is a point,
    # or if WyckoffTranslation lies in the span of TranslationBasis(G).
    # The latter is tested with the Rouche-Capelli theorem.
    if IsSpaceGroup(G) or Length(WyckoffBasis(mywy)) = 0 or
      (RankMatrix(TransposedMat(TranslationBasis(G))) = RankMatrix(TransposedMat(Concatenation(TranslationBasis(G), [WyckoffTranslation(mywy)]))))
      then
      t := WyckoffTranslation(mywy);
    elif Length(TranslationBasis(G)) = 1 or Length(WyckoffBasis(mywy)) = 1 then
      # G is a rod group or Wyckoff position is a line. We have a unique solution.
      # Solve wantedTranslation = currentTranslation + WyckoffBasis*x2 = TranslationBasis * x1
      # (x1|x2) = (TranslationBasis|-WyckoffBasis)^-1 * currentTranslation
      M := TransposedMat(Concatenation(TranslationBasis(G),-WyckoffBasis(mywy)));
      # Since M may be rectangular, we instead compute the pseudoinverse
      # to get the left inverse.
      x := (TransposedMat(M)*M)^-1*TransposedMat(M) * WyckoffTranslation(mywy);
      t := TransposedMat(TranslationBasis(G)) * x{[1..Length(TranslationBasis(G))]};
    else
      # We have a non-unique solution. This could be problematic.
      # The proper thing to do is to find the intersection of several
      # Wyckoff positions (although isn't that just a higher symmetry
      # Wyckoff position?). I'm not sure what a good procedure for that is.
      # However, I don't anticipate creating off-axis groups.
      t := WyckoffTranslation(mywy);
      Error("Situation I can't handle in StandardOriginTransformationOfGroup, a subperiodic group where the WyckoffPosition has an off-axis WyckoffTranslation and no unique choice of way to get it back on-axis.\n");
    fi;
    return TranslationAffineMatrix(-t);
  fi;
  # needtranslate = false, no translation necessary
  return IdentityMat(d+1);
end;

StandardOriginTransformationOfGroupOnScrew := function(G)
  local ODS, op, tloc, P;
  # For L9, L29, L33, which has the origin on a 2-fold screw axis.
  # There are no Wyckoff positions, or at least none that help.
  # So we find the two-fold screw.
  ODS := OrderDetSymMatrixOfGroup(G);
  # Get the first one that has order 2, determinant 1, and is not 
  # symmorphic; get its matrix.
  op := Filtered(ODS, x->x[1]=2 and x[2]=1 and x[3]=false)[1][4];
  # Get the location translation component.
  tloc := DecomposeMatrixTranslationOnRight(TransposedMat(op))[2];
  # Make the translation matrix.
  P := TranslationAffineMatrix(-tloc/2);
  return P;
end;

StandardOriginTransformationOfGroupOnGlide := function(G)
  local ODS, op, tloc, P;
  # pb11 has a sinlge mirror glide plane. Origin is on it.
  # Ditto for p11a
  # So we find the mirror glide.
  ODS := OrderDetSymMatrixOfGroup(G);
  # Get the first one that has order 2, determinant -1, and is not 
  # symmorphic; get its matrix.
  op := Filtered(ODS, x->x[1]=2 and x[2]=-1 and x[3]=false)[1][4];
  # Get the location translation component.
  tloc := DecomposeMatrixTranslationOnRight(TransposedMat(op))[2];
  # Make the translation matrix.
  P := TranslationAffineMatrix(-tloc/2);
  return P;
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
    elif nums = [20] then
      # Another special case. L20 p2_122 have two distinct origin choices. One is conventional
      return StandardOriginTransformationOfGroupOn2z(G);
    elif 18 in nums then
      # L18 has a second lower-symmetry Wyckoff. c2/m11
      nums := FilterLayerGroupNumWithMatchingWyckoff(G, nums);
      if nums = [18] then
        return StandardOriginTransformationOfGroupByWyckoff(G);
      fi;
    elif 24 in nums then
      # L24, pma2, has origin on 2z rather than the mirror.
      nums := FilterLayerGroupNumWithMatchingWyckoff(G, nums);
      if nums = [24] then
        return StandardOriginTransformationOfGroupOn2z(G);
      fi;
    elif nums = [51] or nums = [61] or nums = [63] or nums = [66] or nums = [71,72] or nums = [75] or nums = [80] then
      # p4/m, has a second lower symmetry inversion center.
      # And p4/mbm has inversion centers but one is has a higher order rotation on it.
      # And several trigonal and hexagonal groups have lower symmetry inversion centres.
      return StandardOriginTransformationOfGroupByWyckoff(G);
    elif nums = [52] or nums = [62] or nums = [64] then
      return StandardOriginTransformationOfLayerGroup4FoldTwoChoice(G);
    elif nums = [9] or nums = [33] then
      return StandardOriginTransformationOfGroupOnScrew(G);
    elif 29 in nums then
      if ITNumberOfLayerGroup(G) = 29 then
        return StandardOriginTransformationOfGroupOnScrew(G);
      fi;
    elif nums = [5,12] then
      # p11a and pb11.
      return StandardOriginTransformationOfGroupOnGlide(G);
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
    # (Finding an inversion center is also cheaper than by Wyckoff positions)
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
