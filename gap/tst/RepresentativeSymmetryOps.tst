gap> START_TEST( "RepresentativeSymmetryOps.tst" );
gap> ChangeDirectoryCurrent("..");;
gap> Read("RepresentativeSymmetryOps.gap");;
gap> SetCrystGroupDefaultAction(RightAction);; # Ensure default is respected

# Start with RightAction
gap> TranslationAffineMatrix([1,2,0]);
[[1,0,0,0],[0,1,0,0],[0,0,1,0],[1,2,0,1]]

gap> G := SpaceGroupIT(2,5);; # cm
gap> RepresentativeSymmetryOps(G);
[ [[1,0,0],[0,1,0],[0,0,1]], [[-1,0,0],[0,1,0],[0,0,1]] ]

# Take translations modulo 1 in the representatives to make comparison tractable.
gap> T := [[1,0], [0,1]];;
gap> ops := RepresentativeSymmetryOps(G, T);;
gap> for o in ops do if o[3,1] < 0 then o[3,1] := o[3,1] + 1; fi; if o[3,2] < 0 then o[3,2] := o[3,2] + 1; fi; od;;
gap> ops;
[ [[1,0,0],[0,1,0],[0,0,1]], [[-1,0,0],[0,1,0],[0,0,1]],
  [[1,0,0],[0,1,0],[1/2,1/2,1]], [[-1,0,0],[0,1,0],[1/2,1/2,1]] ]

# Now do LeftAction
gap> SetCrystGroupDefaultAction(LeftAction);;
gap> TranslationAffineMatrix([1,2,0]);
[[1,0,0,1],[0,1,0,2],[0,0,1,0],[0,0,0,1]]

gap> G := SpaceGroupIT(2,5);;
gap> RepresentativeSymmetryOps(G);
[ [[1,0,0],[0,1,0],[0,0,1]], [[-1,0,0],[0,1,0],[0,0,1]] ]

gap> ops := RepresentativeSymmetryOps(G, T);;
gap> for o in ops do if o[1,3] < 0 then o[1,3] := o[1,3] + 1; fi; if o[2,3] < 0 then o[2,3] := o[2,3] + 1; fi; od;;
gap> ops;
[ [[1,0,0],[0,1,0],[0,0,1]], [[-1,0,0],[0,1,0],[0,0,1]],
  [[1,0,1/2],[0,1,1/2],[0,0,1]], [[-1,0,1/2],[0,1,1/2],[0,0,1]] ]

gap> ChangeDirectoryCurrent("tst");
true
gap> STOP_TEST( "RepresentativeSymmetryOps.tst", 10000 );
