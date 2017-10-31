#!/usr/bin/env perl

while(<>) {
    if(/^FLUX /) {
	@line = split();
	$b = $line[1];
	$midpt{$b} += $line[$#line];
	$th{$b} += $line[$#line-1];
	$count{$b}++;
    }
}

foreach $b (sort keys %th) {
    $thave = $th{$b}/$count{$b};
    $midave = $midpt{$b}/$count{$b};
    print "$b $thave $midave\n";
}
