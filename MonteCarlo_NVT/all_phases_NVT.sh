#! /bin/sh

./clean.sh
./Monte_Carlo_NVT.exe "input.gas"
mv output*.0 outputs/NVT/gas
mv config.final outputs/NVT/gas

./Monte_Carlo_NVT.exe "input.liquid"
mv output*.0 outputs/NVT/liquid
mv config.final outputs/NVT/liquid

./Monte_Carlo_NVT.exe "input.solid"
mv output*.0 outputs/NVT/solid
mv config.final outputs/NVT/solid