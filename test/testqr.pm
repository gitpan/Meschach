
my ($ok1, $ok2, $ok3);

my $A= double( pdl [[2,0,1],[0,2,1],[1,0,1]] );
my $b= ~double( pdl [1,4,5] );
my $x=0;

																#### QR Solve : Raw.

my $V ;

my $QR = $A + 0;								# Make a copy

qrfac_($QR,$V);
qrsolve_($x, $b, $QR,$V ); 

my $e = $b - $A x $x ;					# Check result
p "QRsolve Error 1 : $b - $A x $x = $e " unless $ok1= min( abs($e) < 0.000001 );


$ok1
