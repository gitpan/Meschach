
my ($ok1, $ok2, $ok3);

my $A= double( pdl [[2,0,1],[0,2,1],[1,0,1]] );
my $b= ~double( pdl [1,4,5] );
my $x=0;

																#### LU Solve : Raw.

my $Perm ;

my $LU = $A + 0;								# Make a copy

lufac_($LU,$Perm);
lusolve_($x, $b, $LU,$Perm ); 

my $e = $b - $A x $x ;					# Check result
p "LUsolve Error 1 : $b - $A x $x = $e " unless $ok1= min( $e == 0 );


																#### LU Solve : Friendly.

($LU,$Perm,$x) = lusolve( $b, $A );
$e = $b - $A x $x ;

p "LUsolve Error 2 : $b - $A x $x = $e " unless $ok2= min( $e == 0 );



$x = lusolve( $b, $LU, $Perm );
$e = $b - $A x $x ;

p "LUsolve Error 3 : $b - $A x $x = $e " unless $ok3= min( $e == 0 );

$ok1 && $ok2 && $ok3; 

