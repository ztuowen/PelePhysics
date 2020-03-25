
CHEMINP=chem.inp
THERMINP=therm.dat
FINALFILE=jp8-lu-red.cpp

FMC=${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/products/bin/fmc.py
HEADERDIR=${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/header

echo ${FUEGO_PYTHON}

${FUEGO_PYTHON} ${FMC} -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${FINALFILE}

echo Compiling ${FINALFILE}...
