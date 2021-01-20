#!/bin/sh

./clean.sh
./MolecularDynamics "input.gas"
mv output*.dat gas_phase
mv ave_*.out gas_phase
mv config.final gas_phase

./MolecularDynamics "input.liquid"
mv output*.dat liquid_phase
mv ave_*.out liquid_phase
mv config.final liquid_phase

./MolecularDynamics "input.solid"
mv output*.dat solid_phase
mv ave_*.out solid_phase
mv config.final solid_phase