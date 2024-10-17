alltests := DirectoryCurrent();

TestDirectory( alltests,
  rec( exitGAP := true, exclude := [],
       testOptions := rec( compareFunction := "uptowhitespace") ) );

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
