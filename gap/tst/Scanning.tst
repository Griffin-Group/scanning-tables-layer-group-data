gap> START_TEST("Scanning.tst");
gap> ChangeDirectoryCurrent("..");;
gap> Read("Scanning.gap");;

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
