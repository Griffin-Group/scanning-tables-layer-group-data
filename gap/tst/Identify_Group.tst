gap> START_TEST( "Identify_Group.tst" );
gap> ChangeDirectoryCurrent("..");;
gap> Read("Identify_Group.gap");;
gap> SetCrystGroupDefaultAction(RightAction);;

# Step 1, check that base functions work
gap> OrderDetOfMatrix([[1,0],[0,1]]);
[1,1]
gap> OrderDetOfMatrix([[0,1],[-1,0]]);
[4,1]
gap> OrderDetOfMatrix([[-1,0],[0,-1]]);
[-1,-1]
gap> OrderDetOfMatrix([[-1,0],[0,1]]);
[2,-1]

gap> OrderDetSymOfMatrix(IdentityMat(3));
[1,1,true]
gap> OrderDetSymOfMatrix([[-1,0,0],[0,1,0],[1/2,0,1]]);
[2,-1,true]
gap> OrderDetSymOfMatrix([[-1,0,0],[0,1,0],[0,1/2,1]]);
[2,-1,false]

gap> OrderDetSymMatrix(IdentityMat(3));
[1,1,true,[[1,0,0],[0,1,0],[0,0,1]]]

gap> DecomposeMatrixTranslationOnRight([[-1,0,0],[0,1,0],[1/2,0,1]]);
[[0,0],[1/2,0]]
gap> DecomposeMatrixTranslationOnRight([[-1,0,0],[0,1,0],[0,1/2,1]]);
[[0,1/2],[0,0]]

gap> G := SpaceGroupIT(2,3);; # pm
gap> IsMatrixSymmorphicEquivalent(IdentityMat(3), G);
true
gap> IsMatrixSymmorphicEquivalent([[-1,0,0],[0,1,0],[0,1,1]], G);
true

gap> OrderDetMatrixOfGroup(G);
[ [1,1,[[1,0,0],[0,1,0],[0,0,1]]], [2,-1,[[-1,0,0],[0,1,0],[0,0,1]]] ]
gap> OrderDetSymMatrixOfGroup(G);
[ [1,1,true,[[1,0,0],[0,1,0],[0,0,1]]], [2,-1,true,[[-1,0,0],[0,1,0],[0,0,1]]] ]
gap> OrderDetSymOfGroup(G);
[ [1,1,true], [2,-1,true] ]

gap> G := SpaceGroupIT(2,4);; # pg
gap> IsMatrixSymmorphicEquivalent([[-1,0,0],[0,1,0],[0,1/2,1]], G);
false
gap> OrderDetSymOfGroup(G);
[ [1,1,true], [2,-1,false] ]

gap> CrossProduct([1,0,0], [0,1,0]);
[0,0,1]

# Now do the real test. Identify group numbers
gap> List([1..75], i -> ITNumberOfRodGroup(RodGroups[i])) = [1..75];
true

gap> C := [[1, -2, 3, 0], [1, 0, -1, 0], [-2, 3, 1, 0], [1/2,-2,3/5,1]];;
gap> List([1..75], i -> ITNumberOfRodGroup(RodGroups[i]^C)) = [1..75];
true

gap> List([1..80], i -> ITNumberOfLayerGroup(LayerGroups[i])) = [1..80];
true

gap> List([1..80], i -> ITNumberOfLayerGroup(LayerGroups[i]^C)) = [1..80];
true

# Identify settings
# NB. Not all orthorhombic rod groups have all 6 settings distinct.
gap> for i in [1,2,3,4,5,6,7,17,18,19] do for setting in SubPeriodicGroupSettingsIT("Rod", i) do if (SettingOfRodGroupNC(RodGroupIT(i, setting), i) <> setting) then Print("Rod ", i, " setting ", setting, " falsely ID'd.\n"); fi; od; od;;

gap> for i in Concatenation([8..16],[20..22]) do for setting in ["abc","cba","bca"] do if (SettingOfRodGroupNC(RodGroupIT(i, setting), i) <> setting) then Print("Rod ", i, " setting ", setting, " falsely ID'd.\n"); fi; od; od;;

# Layer groups
gap> for i in [1..80] do for setting in SubPeriodicGroupSettingsIT("Layer", i) do if (SettingOfLayerGroupNC(LayerGroupIT(i, setting), i) <> setting) then Print("Layer ", i, " setting ", setting, " falsely ID'd.\n"); fi; od; od;;

# Test the other group action
gap> SetCrystGroupDefaultAction(LeftAction);;

gap> OrderDetSymOfMatrix(IdentityMat(3));
[1,1,true]
gap> OrderDetSymOfMatrix([[-1,0,1/2],[0,1,0],[0,0,1]]);
[2,-1,true]
gap> OrderDetSymOfMatrix([[-1,0,0],[0,1,1/2],[0,0,1]]);
[2,-1,false]

gap> G := SpaceGroupIT(2,4);; # pg
gap> IsMatrixSymmorphicEquivalent([[-1,0,0],[0,1,1/2],[0,0,1]], G);
false
gap> OrderDetSymOfGroup(G);
[ [1,1,true], [2,-1,false] ]

gap> List([4,14,24,34,44,54,64,74], i -> ITNumberOfRodGroup(RodGroupIT(i)));
[4,14,24,34,44,54,64,74]
gap> List([4,14,24,34,44,54,64,74], i -> ITNumberOfRodGroup(RodGroupIT(i)^TransposedMat(C)));
[4,14,24,34,44,54,64,74]
gap> List([4,14,24,34,44,54,64,74], i -> ITNumberOfLayerGroup(LayerGroupIT(i)));
[4,14,24,34,44,54,64,74]
gap> List([4,14,24,34,44,54,64,74], i -> ITNumberOfLayerGroup(LayerGroupIT(i)^TransposedMat(C)));
[4,14,24,34,44,54,64,74]

gap> SettingOfRodGroupNC(RodGroupIT(6, "bac"),6);
"bac"
gap> SettingOfLayerGroupNC(LayerGroupIT(24,'b'),24);
'b'

# A tricky case that came up.
gap> G := LayerGroupIT(36) ^ [[-1,-1,0,0],[1,-1,0,3/4],[0,0,1,0],[0,0,0,1]];;
gap> ITNumberOfLayerGroup(G);
36

gap> ChangeDirectoryCurrent("tst");
true
gap> STOP_TEST( "Identify_Group.tst" );
