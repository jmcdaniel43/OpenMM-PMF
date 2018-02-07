#!/usr/bin/perl
#
# packmol doesn't put box vectors into pdb input file, so that's what this script does
# add 2-3 angstrom buffer to packmol box.  packmol doesn't consider pbc, so if we don't
# add a buffer there could be very close/high energy contacts accross pbc.

$buffer=3.0;
$pdbfile=$ARGV[0];
$packmolinput=$ARGV[1];

$tempfile="temppp";
# get box vector from packmolinput

open ifile, $packmolinput, or die;

while(<ifile>)
{ 
if (/inside cube/){last}
}
@array=split;$box=$array[5];
$box=$box+$buffer;
close(ifile);

# now open pdb and print box
open ifile, $pdbfile, or die;
open (ofile,">$tempfile");

$flag=0;
while(<ifile>)
{
if (/HETATM/){
if ( $flag == 0 ){
#print box
printf ofile ("%6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s\n","CRYST1", $box, $box, $box, 90, 90, 90," P 1           1");
$flag=1;
}
}

print ofile;
}

rename( $tempfile ,  $pdbfile);
