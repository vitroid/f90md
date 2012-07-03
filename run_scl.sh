#!/bin/sh

cat <<EOF > test.input
[INTRPAIR]
1 1 LJ
0.99768d0 3.41d0
[STEPS]
100000
[INTVPS]
0.001d0
[LOGINTV]
1000
EOF
./scl >> test.input
./main < test.input

