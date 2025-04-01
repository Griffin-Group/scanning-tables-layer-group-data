gap> START_TEST("Scanning.tst");
gap> ChangeDirectoryCurrent("..");;
gap> Read("Scanning.gap");;

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

# Test sectional group methods
# p112 -> p121 or p1
gap> G := SectionalGroupOfLayerGroup(LayerGroupIT(3), [1,0,0], [0,1,0], 0);
<matrix group with 3 generators>
gap> GeneratorsOfGroup(G);
[ [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] ], 
  [ [ -1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, -1, 0 ], [ 0, 0, 0, 1 ] ], 
  [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 1 ], [ 0, 0, 0, 1 ] ] ]
gap> G := SectionalGroupOfLayerGroup(LayerGroupIT(3), [1,0,0], [0,1,0], 1/4);
<matrix group with 2 generators>
gap> GeneratorsOfGroup(G);
[ [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] ], 
  [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 1 ], [ 0, 0, 0, 1 ] ] ]

gap> ScanningVectorFromIntVector([1,0,0]);
[0,1,0]
gap> ScanningVectorFromIntVector([0,1,0]);
[-1,0,0]
gap> ScanningVectorFromIntVector([1,-1,0]);
[1,0,0]

gap> G := LayerGroupIT(58);;
gap> ScanningVectorFromIntVectorAndGroup([1,0,0], G);
[0,1,0]
gap> ScanningVectorFromIntVectorAndGroup([1,-1,0], G);
[1,1,0]

gap> G := LayerGroupIT(35);; # cm2m
gap> ScanningVectorFromIntVectorAndGroup([1,0,0], G);
[0,1,0]
gap> d := ScanningVectorFromIntVectorAndGroup([2,1,0], G);;
gap> 2 * d[2] - 1 * d[1]; # Cross product
1/2
gap> d[3];
0

gap> G := LayerGroupIT(27);;
gap> H := ScanningGroupOfLayerGroup(G, [1,0,0]);;
gap> H = G;
true
gap> ScanningGroupOfLayerGroup(G, [0,1,0]) = LayerGroupIT(27, 'b');
true
gap> ScanningGroupOfLayerGroup(G, [1,1,0]) = LayerGroupIT(4);
true
gap> G := LayerGroupIT(7);;
gap> ScanningGroupOfLayerGroup(G, [1,0,0]) = G;
true
gap> ScanningGroupOfLayerGroup(G, [0,1,0]) = LayerGroupIT(7, '3');
true
gap> ScanningGroupOfLayerGroup(G, [0,1,0], [-1,-1,0]) = LayerGroupIT(7, '2');
true

gap> SpecialScanSOfLayerGroup(LayerGroupIT(9), [1,0,0], [0,1,0]);
[0, 1/2]
gap> SpecialScanSOfLayerGroup(LayerGroupIT(9), [0,1,0], [-1,0,0]);
[ ]

gap> data := ScanningTableFun(LayerGroupIT(9), [1,0,0], [0,1,0]);;
gap> reference := rec(H := 9, H_origin:=[0,0,0], H_setting:='a',
>      d:=[0,1,0], general:=1, general_length:=1, general_origin:=[0,0,0], general_setting:='1',
>      special := [rec(length:=1, number:=9, origin:=[0,0,0], s:=[0,1/2], setting:="abc")]);;
gap> data = reference;
true
gap> ScanningTable(LayerGroupIT(9), [1,0,0]) = reference;
true
gap> ScanningTableFromLayerIT(9, [1,0,0]) = reference;
true

gap> reference.orientation := [[1,0,0]];;
gap> reference.d := [[0,1,0]];;
gap> reference2 := rec(H:=9, H_origin:=[0,0,0], H_setting:='b', d:=[[-1,0,0]], orientation:=[[0,1,0]],
>     general:=1, general_length:=1, general_origin:=[0,0,0], general_setting:='1', special:=[]);;
gap> FullScanningTableFromLayerIT(9) = [reference, reference2];
Scanning layer # 9 and orientation [ 1, 0, 0 ].
Scanning layer # 9 and orientation [ 0, 1, 0 ].
true

gap> FindLinearOrbitsFromLayerIT(5);
rec( a := [ [ 0 ], [ 1/4 ], [ 1/2 ], [ 3/4 ], [ 1/7 ] ], 
  b := [ [ 0, 1/2 ], [ 1/4, 3/4 ], [ 1/7, 9/14 ] ] )

gap> FindLinearOrbitsFromLayerIT(10);
rec( a := [ [ 0, 1/2 ], [ 1/4, 3/4 ], [ 1/7, 5/14, 9/14, 6/7 ] ], 
  ao := [ [ 0, 1/2 ], [ 1/4, 3/4 ], [ 1/7, 5/14, 9/14, 6/7 ] ], 
  b := [ [ 0, 1/2 ], [ 1/4, 3/4 ], [ 1/7, 9/14 ] ], 
  bo := [ [ 0, 1/2 ], [ 1/4, 3/4 ], [ 1/7, 9/14 ] ] )

gap> ChangeDirectoryCurrent("tst");
true
gap> STOP_TEST("Scanning.tst");
