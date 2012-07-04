#!/bin/sh

#100+100steps
cat <<EOF > @2compo.1.input
[INTRPAIR]
1 1 LJ
0.99768d0 3.41d0
[INTRPAIR]
2 2 LJ
0.99768d0 3.41d0
[INTRPAIR]
1 2 LB
[STEPS]
100
[INTVPS]
0.001d0
[LOGINTV]
10
EOF
./gen_2compo >> @2compo.1.input
./main < @2compo.1.input > @2compo.1.out
./main < @2compo.1.out > @2compo.2.out

#200steps at once
cat <<EOF > @2compo.1+2.input
[INTRPAIR]
1 1 LJ
0.99768d0 3.41d0
[INTRPAIR]
2 2 LJ
0.99768d0 3.41d0
[INTRPAIR]
1 2 LB
[STEPS]
200
[INTVPS]
0.001d0
[LOGINTV]
10
EOF
./gen_2compo >> @2compo.1+2.input
./main < @2compo.1+2.input > @2compo.1+2.out
