#!/usr/local/bin/perl

#Angstroem to nm
$DIV = 10;
 
$xxx = @ARGV;
if ($xxx != 3) {
                die "usage: noe_inp.pl <pdb_noe.file> <g96 cor file> <output file for unparsed NOEs>\n";
               }
$INNOE  = $ARGV[0];
$INTOPO = $ARGV[1];
$REST   = $ARGV[2];
 
open(DATA, $INNOE);
 
$num=1;
while ( <DATA> )
{
 $line= $_;
 chomp($_);
 @line = split (/\s+/,$_);
                       $a[$num] = $line[1];
                       $XX[$num]= $line[1]; 
                       $b[$num] = $line[2];
                       $c[$num] = $line[3];
                       $d[$num] = $line[4];
                       $e[$num] = $line[5];
                       $f[$num] = $line[6];
                       $g[$num] = $line[7];
                       $num++;
}
close (DATA);

#Make Substitutions on noe list
@RES = (ASN, CYS, HIS,HIS,ILE,"LYS+",PHE,PHE,TYR,TYR);
@ORI = (QD2,  HG, HE1,HD2,QD1, HE3  , QD, QE, QD, QE);
@SUB = (AQD2, CH,HHE1,HD1,IQD, LH3  ,PQD,PQE,TQD,TQE);

for ($y=1; $y < @a ; $y++)
 {
  for ($l=0; $l < @RES; $l++)
   {
    if ($b[$y] eq $RES[$l] && $c[$y] eq $ORI[$l])
     {
      $c[$y] = $SUB[$l];
     }
    if ($e[$y] eq $RES[$l] && $f[$y] eq $ORI[$l])
     {
      $f[$y] = $SUB[$l];
     }
    else {next}
   }
 }

open(TOP, $INTOPO);
$numt=1;
while ( <TOP> )
{
 $line= $_;
 @line = split (/\s+/,$_);
                       $aa[$numt] = $line[1];
                       $bb[$numt] = $line[2];
                       $cc[$numt] = $line[3];
                       $dd[$numt] = $line[4];
                       $numt++;
 }
close (TOP);

print "TITLE \n";
print "NOE specification file for $INNOE \n";
print "END \n";
print "NOE \n";           

open(OUT, ">$REST"); 

#NOE DEFINITION LIST

@NOE = (HN,HA,QB,HA1,HA2,QQD,CH,HG,HB,HD1,HE1,HZ2,HH2,HZ3,QG1,QG2,QG,QD,QA,HZ,QQG,HE3,HB2,HB3,QD1,QD2,IQD,QE,HG2,HG3,HG12,HG13,PQD,TQD,PQE,TQE,HD21,HD22,HD2,         HD3,AQD2,HE21,HE22, QE2,HHE1,HE2,LH3);
@TOP = ( H,CA,CB, CA, CA,CD1,HG,CG,CB,HD1,HE1,HZ2,HH2,HZ3,CG1,CG2,CG,CD,CA,HZ,CG1,HE3, CB, CB,CD1,CD2, CD,CE, CG, CG, CG1, CG1, CG, CG, CZ, CZ,HD21,HD22, CD,          CD,HD21,HE21,HE22,HE21, CE1, CE, CE);
@DEF = ( E, S, N,  N,  N,  S, E, N, E,  E,  E,  E,  E,  E,  N,  N, N, N, N, E,  S,  E,  S,  S,  N,  N,  N, N,  S,  S,   S,   S,  A,  A,  A,  A,   E,   E,  N,          N,    N,   E,   E,   N,   E,  N,  N);
@STE = ( 0, 0, 0,  0,  0,  1, 0, 0, 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0,  1,  0,  0,  0,  0,  0,  0, 0,  0,  0,   0,   0,  0,  0,  0,  0,   0,   0,  0,          0,    1,   0,   0,   1,   0,  0,  0);

#Start Parsing

for ($i=1; $i < @c; $i++)
 {
  for ($q=0; $q < @NOE; $q++)
   {
    for ($r=0; $r < @NOE; $r++)
     {
      if ($c[$i] eq $NOE[$q] && $f[$i] eq $NOE[$r])
       {
        if ($STE[$q] != 0 && $STE[$r] == 0)
        {
         $A1=$TOP[$q];
         $A2=$TOP[$r];
         $TYPE1=$DEF[$q];
         $TYPE2=$DEF[$r];
         $XX[$i] = X;
         &PARSEV 
        }
        if ($STE[$r] != 0 && $STE[$q] == 0)
        {
         $A1=$TOP[$q];
         $A2=$TOP[$r];
         $TYPE1=$DEF[$q];
         $TYPE2=$DEF[$r];
         $XX[$i] = X;
         &PARSEV 
        }
        if ($STE[$r] != 0 && $STE[$q] != 0)
        {
         $A1=$TOP[$q];
         $A2=$TOP[$r];
         $TYPE1=$DEF[$q];
         $TYPE2=$DEF[$r];
         $XX[$i] = X;
         &PARSEV
        }
        if ($STE[$q] == 0 && $STE[$r] == 0) 
        {
         $A1=$TOP[$q];
         $A2=$TOP[$r];
         $TYPE1=$DEF[$q];
         $TYPE2=$DEF[$r];
         $XX[$i] = X;
         &PARSE
        }                                                
       }
   #   else 
   #    {
   #     next;
   #    }
     }
   }
 }
sub PARSE
{
 for ($f=1; $f < @aa; $f++)
  {
   if ($aa[$f] eq $a[$i] && $cc[$f] eq $A1)
    {
     print "$dd[$f]$TYPE1  ";
     for ($g=1; $g < @aa; $g++)
      {
       if ($aa[$g] eq $d[$i] && $cc[$g] eq $A2)
        {
         $nm=$g[$i]/$DIV;
         print "$dd[$g]$TYPE2  $nm\n";
        }
      }
    }
  }
}                             

sub PARSEV
{
 if ($STE[$q] != 0 && $STE[$r] == 0)
  {
   for ($f=1; $f < @aa; $f++)
    {
     if ($aa[$f] eq $a[$i] && $cc[$f] eq $A1)
      {
       print "$dd[$f]$TYPE1$dd[$f+$STE[$q]]  ";
       for ($g=1; $g < @aa; $g++)
        {
         if ($aa[$g] eq $d[$i] && $cc[$g] eq $A2)
          {
           $nm=$g[$i]/$DIV;
           print "$dd[$g]$TYPE2  $nm\n";
          }
        }
      }
    }
  }
 if ($STE[$r] != 0 && $STE[$q] == 0)
  { 
   for ($f=1; $f < @aa; $f++)
    {
     if ($aa[$f] eq $a[$i] && $cc[$f] eq $A1)
      {
       print "$dd[$f]$TYPE1  ";
       for ($g=1; $g < @aa; $g++)
        {
         if ($aa[$g] eq $d[$i] && $cc[$g] eq $A2)
          {
           $nm=$g[$i]/$DIV;
           print "$dd[$g]$TYPE2$dd[$g+$STE[$r]]  $nm\n";
          }
        }
      }
    }
  }
 if ($STE[$r] != 0 && $STE[$q] != 0)
  {
   for ($f=1; $f < @aa; $f++)
    {
     if ($aa[$f] eq $a[$i] && $cc[$f] eq $A1)
      {
       print "$dd[$f]$TYPE1$dd[$f+$STE[$q]]  ";
       for ($g=1; $g < @aa; $g++)
        {
         if ($aa[$g] eq $d[$i] && $cc[$g] eq $A2)
          {
           $nm=$g[$i]/$DIV;
           print "$dd[$g]$TYPE2$dd[$g+$STE[$r]]  $nm\n";
          }
        }
      }
    }
  }
}

for ($s=1; $s < @a ; $s++)
 {
  if ($XX[$s] ne X)
   {
     print OUT "$XX[$s] $b[$s] $c[$s] $d[$s] $e[$s] $f[$s] \n";
   }
 }
close (OUT);
print "END \n";                   
