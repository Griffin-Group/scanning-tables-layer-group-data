#!/usr/bin/env perl
use strict;
use warnings;

# Open the template
open(my $in, "<", "subgrp_template.grp") or die;
read($in, my $template, -s $in);
close $in or die;
# Go over the different files
foreach ("frieze.gi","rod.gi","layer.gi") {
	# Open (slurp) the file
	open(my $in, "<", $_) or die;
	read($in, my $contents, -s $in);
	close $in or die;
	# Clean up some of the whitespace
	$contents =~ s/([,\[]) /$1/g;
	$contents =~ s/ ([\]\)])/$1/g;
	$contents =~ s/^\s+//gm;
	# Insert into the template
	$template =~ s/$_/$contents/;
}
# Output
open(my $out, ">", "subgrp.grp") or die;
print $out $template;
close $out or die;
