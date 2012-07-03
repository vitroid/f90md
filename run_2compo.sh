#!/bin/sh

cat <<EOF > 2compo.input
[INTRPAIR]
1 1 LJ
0.99768d0 3.41d0
[INTRPAIR]
2 2 LJ
0.99768d0 3.41d0
[INTRPAIR]
1 2 LB
[STEPS]
100000
[INTVPS]
0.001d0
[LOGINTV]
1000
EOF
./gen_2compo >> 2compo.input
./main < 2compo.input

