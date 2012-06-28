#!/bin/sh

cp test.head test.input
./scl >> test.input
./main < test.input

