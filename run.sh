#!/bin/sh

seeds=(7 22 100 222 273687)

for seed in ${seeds[*]}
do
    ./bin/practica1 zoo_set.dat zoo_set_const_10.const 7 $seed
    ./bin/practica1 zoo_set.dat zoo_set_const_20.const 7 $seed
    ./bin/practica1 glass_set.dat glass_set_const_10.const 7 $seed
    ./bin/practica1 glass_set.dat glass_set_const_20.const 7 $seed
    ./bin/practica1 bupa_set.dat bupa_set_const_10.const 16 $seed
    ./bin/practica1 bupa_set.dat bupa_set_const_20.const 16 $seed
done