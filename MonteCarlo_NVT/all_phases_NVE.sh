#! /bin/sh

./MolDyn/clean.sh
./MolDyn/MolecularDynamics "MolDyn/input.gas"
mv ave*.out outputs/NVE/gas

./MolDyn/clean.sh
./MolDyn/MolecularDynamics "MolDyn/input.liquid"
mv ave*.out outputs/NVE/liquid

./MolDyn/clean.sh
./MolDyn/MolecularDynamics "MolDyn/input.solid"
mv ave*.out outputs/NVE/solid