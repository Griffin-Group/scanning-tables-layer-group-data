#!/bin/env gap


TranslationAffineMatrix := function(v)
  local d, op;
  #% Creates an affine matrix with a translation v
  if not IsVector(v) then;
    Error("Expected v to be a vector.");
  fi;
  # Get the dimensionality
  d := Length(v);
  op := IdentityMat(d+1);
  if CrystGroupDefaultAction = LeftAction then;
    op{[1..d]}[d+1] := v;
  else;
    op[d+1]{[1..d]} := v;
  fi;
  return op;
end;

# Let us collect the representative symmetry operations for p6/mmm
# This comes from a coset decomposition with respect to the translation subgroup
TranslationCrystGroupByBasis := function(B)
  local ops, d;
  #% Turns a basis B (List of vectors) into an AffineCrystGroup
  if not IsMatrix(B) then
    Error("Expected B to be a matrix.");
  fi;
  # Get dimension
  d := Length(B[1]);
  ops := List(B, TranslationAffineMatrix);
  return AffineCrystGroup(ops);
end;

TranslationSubgroup := function(G)
  local T;
  #% Gets the translation subgroup of AffineCrystGroup G
  if not (IsAffineCrystGroupOnLeft(G) or IsAffineCrystGroupOnRight(G)) then
    Error("G not an AffineCrystGroup");
  fi;
  T := TranslationCrystGroupByBasis(TranslationBasis(G));
  if (IsAffineCrystGroupOnLeft(G) and CrystGroupDefaultAction = RightAction) or
    (IsAffineCrystGroupOnRight(G) and CrystGroupDefaultAction = LeftAction) then
    T := TransposedMatrixGroup(T);
  fi;
  return T;
end;


RepresentativeSymmetryOps := function(G, T...)
  local translation, cosets, ops, op, g;
  #% Get the representative symmetry operations of some crystal group G,
  #% relative to its translation basis. i.e. the coset decomposition of G
  #% with respect to T.
  #% You can provide a translation subgroup T if you want to consider a
  #% non-primitive unit cell. You may pass a list of vectors (i.e. a basis)
  #% or a Group.
  if not (IsAffineCrystGroupOnLeft(G) or IsAffineCrystGroupOnRight(G)) then
    Error("Expected G to be an AffineCrystGroup.");
  fi;
  # Extract the main translation group
  if Length(T) > 0 then;
    if IsGroup(T[1]) and IsAffineCrystGroup(T[1]) then;
      translation := T[1];
    else;
      # T can be a list of basis vectors too
      translation := TranslationCrystGroupByBasis(T[1]);
    fi;
    # Verify that this is indeed a translation subgroup of G
    if not IsSubset(G, translation) then;
      Error("T is not a subgroup of G.");
    fi;
  else;
    translation := TranslationSubgroup(G);
  fi;
  # Get the coset decomposition
  cosets := LeftCosets(G, translation);
  # Get the Representatives from this list
  ops := [];
  for g in cosets do
    op := MutableCopyMat(Representative(g));
    Add(ops, op);
  od;
  return ops;
end;
