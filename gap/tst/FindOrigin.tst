gap> START_TEST("FindOrigin.tst");
gap> ChangeDirectoryCurrent("..");;
gap> Read("FindOrigin.gap");;

gap> SetCrystGroupDefaultAction(LeftAction);;

# Test origin identification
# Rod group with inversion
gap> G := RodGroupIT(12);;
gap> T := [[1,0,0,1/2],[0,1,0,-1/2],[0,0,1,1/4],[0,0,0,1]];;
gap> StandardOriginTransformationOfGroup(G) = IdentityMat(4);
true
gap> StandardOriginTransformationOfGroup(G^T) = Inverse(T);
true
gap> OriginLocationOfGroup(G);
[0,0,0]
gap> OriginLocationOfGroup(G^T);
[1/2,-1/2,1/4]
gap> OriginLocationOfGroup(G^[[1,0,0,0],[0,1,0,0],[0,0,1,1/2],[0,0,0,1]]);
[0,0,0]

# Rod group without inversion
gap> G := RodGroupIT(13);;
gap> T := [[1,0,0,1/2],[0,1,0,-1/2],[0,0,1,1/4],[0,0,0,1]];;
gap> StandardOriginTransformationOfGroup(G) = IdentityMat(4);
true
gap> StandardOriginTransformationOfGroup(G^T) = Inverse(T);
true
gap> OriginLocationOfGroup(G);
[0,0,0]
gap> OriginLocationOfGroup(G^T);
[1/2,-1/2,1/4]

# Special case: R14 and R19
gap> OriginLocationOfGroup(RodGroupIT(19));
[0,0,0]
gap> OriginLocationOfGroup(RodGroupIT(19)^[[1,0,0,0],[0,1,0,0],[0,0,1,1/4],[0,0,0,1]]);
[0,0,1/4]
gap> OriginLocationOfGroup(RodGroupIT(19)^[[1,0,0,0],[0,1,0,0],[0,0,1,3/4],[0,0,0,1]]);
[0,0,1/4]
gap> OriginLocationOfGroup(RodGroupIT(14));
[0,0,0]
gap> OriginLocationOfGroup(RodGroupIT(14)^[[1,0,0,0],[0,1,0,0],[0,0,1,1/4],[0,0,0,1]]);
[0,0,1/4]

# Special case: L47 (one inversion center is not the origin)
gap> G := LayerGroupIT(47);;
gap> T := [[1,0,0,1/4],[0,1,0,1/4],[0,0,1,0],[0,0,0,1]];;
gap> StandardOriginTransformationOfGroup(G) = IdentityMat(4);
true
gap> StandardOriginTransformationOfGroup(G^T) = Inverse(T);
true
gap> OriginLocationOfGroup(G);
[0,0,0]
gap> OriginLocationOfGroup(G^T);
[1/4,1/4,0]

# Special case: L48 (setting-dependent origin)
gap> G := LayerGroupIT(48);;
gap> OriginLocationOfGroup(G);
[0,0,0]
gap> OriginLocationOfGroup(G^[[0,1,0,0],[-1,0,0,0],[0,0,1,0],[0,0,0,1]]);
[1/4,1/4,0]

# Special case: tetragonal groups with two origin choices.
gap> OriginLocationOfGroup(LayerGroupIT(62,'1'));
[1/4,1/4,0]
gap> OriginLocationOfGroup(LayerGroupIT(62,'2'));
[0,0,0]

# Special case: L20 p2_122
# (If you don't catch this case properly, you'll get [0,0,0] when shifted.)
gap> G := LayerGroupIT(20);;
gap> OriginLocationOfGroup(G);
[0,0,0]
gap> OriginLocationOfGroup(G^[[1,0,0,1/4],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/4,0,0]
gap> OriginLocationOfGroup(G^[[0,-1,0,3/8],[-1,0,0,3/4],[0,0,-1,0],[0,0,0,1]]);
[3/8,3/4,0]

# Special case: L18 c2/m11.
# It has a lower-symmetry inversion center.
gap> G := LayerGroupIT(18);;
gap> OriginLocationOfGroup(G^[[1,0,0,1/4],[0,1,0,1/4],[0,0,1,0],[0,0,0,1]]);
[1/4,1/4,0]

# Case: L24, pma2. Origin is on 2, not m.
gap> G := LayerGroupIT(24);;
gap> OriginLocationOfGroup(G^[[1,0,0,1/4],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/4,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(31)^[[1,0,0,1/4],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/4,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(36)^[[1,0,0,1/4],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/4,1/4,0]

# Some cases where the origin lies on a screw or glide with no other reference.
gap> OriginLocationOfGroup(LayerGroupIT(9)^[[1,0,0,0],[0,1,0,1/4],[0,0,1,0],[0,0,0,1]]);
[0,1/4,0]
gap> OriginLocationOfGroup(LayerGroupIT(9)^[[1,0,0,0],[0,1,0,1/3],[0,0,1,0],[0,0,0,1]]);
[0,1/3,0]
gap> OriginLocationOfGroup(LayerGroupIT(12)^[[1,0,0,1/4],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/4,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(29)^[[1,0,0,1/4],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/4,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(33)^[[1,0,0,1/4],[0,1,0,0],[0,0,1,0],[0,0,0,1]]); # This one also has the origin on the screw but not glide
[1/4,0,0]

# Check the 4-fold behaves as expected
gap> OriginLocationOfGroup(LayerGroupIT(49)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]

# Ah, here we have another case of multiple inversion centers, some of lower symmetry.
gap> OriginLocationOfGroup(LayerGroupIT(51)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]

# L54, p42_12 has two similarly-weighted Wyckoffs. Make sure we get the 4-fold axis.
# I know we start on the 4-fold axis. There's a few others with the same issue.
gap> G := LayerGroupIT(54);;
gap> OriginLocationOfGroup(G);
[0,0,0]
gap> OriginLocationOfGroup(G^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(56)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(58)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(60)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(63)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]

# L64, has a subtle setting switch between origins in origin choice 2.
gap> OriginLocationOfGroup(LayerGroupIT(64));
[0,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(64)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(64)^[[1,0,0,1/2],[0,1,0,1/2],[0,0,1,0],[0,0,0,1]]);
[0,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(64)^[[0,1,0,0],[-1,0,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(64)^[[1/2,1/2,0,0],[-1/2,1/2,0,0],[0,0,1,0],[0,0,0,1]]);
[0,0,0]

# L62, has a subtle setting switch between origins in origin choice 2.
gap> OriginLocationOfGroup(LayerGroupIT(62));
[0,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(62)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(62)^[[1,0,0,1/2],[0,1,0,1/2],[0,0,1,0],[0,0,0,1]]);
[0,0,0]

# L52, has a subtle setting switch between origins in origin choice 2.
gap> OriginLocationOfGroup(LayerGroupIT(52));
[0,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(52)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(52)^[[1,0,0,1/2],[0,1,0,1/2],[0,0,1,0],[0,0,0,1]]);
[0,0,0]

# The following groups have a lower symmetry inversion center which we can get mixed up with.
gap> OriginLocationOfGroup(LayerGroupIT(66)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(71)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(72)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(75)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]
gap> OriginLocationOfGroup(LayerGroupIT(80)^[[1,0,0,1/2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
[1/2,0,0]

gap> ChangeDirectoryCurrent("tst");
true
gap> STOP_TEST("FindOrigin.tst");