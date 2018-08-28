#!/bin/bash

awk ' { t = $2; $2 = $14; $14 = t; print; } ' $1 > temp.txt

awk ' { t = $3; $3 = $16; $16 = t; print; } ' temp.txt > temp2.txt
awk ' { t = $8; $8 = $11; $11 = t; print; } ' temp2.txt > $2 
