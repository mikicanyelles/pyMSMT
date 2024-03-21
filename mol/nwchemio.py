"""
Module for writting a NWChem file and read the coordinates and force
constants from NWChem output file.
"""
from __future__ import absolute_import, print_function, division
import numpy
import linecache
#from pymsmtexp import *
#from msmtmol.constants import B_TO_A
#from msmtmol.pt import AtomicNum

from exp import *
from mol.constants import B_TO_A
#from mol.periodic_table import AtomicNum

#------------------------------------------------------------------------------
#--------------------------Write GAMESS input file-----------------------------
#------------------------------------------------------------------------------


def write_nwchem_optf(nwcoptf2, smchg, SpinNum, gatms, signum=3, max_iter=1000):
    """
    Writes an input file for a geometry optimization calculation in NWChem.

    Parameters:
        goptf2 (str): The name of the input file to write.
        smchg (int): The charge of the system.
        SpinNum (int): The multiplicity of the system.
        gatms (list): A list of `MoleculeAtom` objects representing the atoms in the system.
        signum (int): The number of decimal places to use for the atomic coordinates.
        max_iter (int): The maximum number of iterations for the optimization.
    """
    with open(nwcoptf2, 'w') as f:
        f.write('start small_model\n')
        f.write('title "Geometry Optimization"\n')
        f.write('charge {}\n'.format(smchg))
        f.write('basis\n')
        f.write('  * library 6-31g*\n')
        f.write('end\n')
        f.write('geometry\n')
        for gatm in gatms:
            element = gatm.element
            #nuchg = AtomicNum[element]
            #nuchg = round(nuchg, signum)
            xyz = gatm.crdx, gatm.crdy, gatm.crdz
            if signum == 3:
                f.write("  {:<2} {:10.3f}{:10.3f}{:10.3f}\n".format(element, *xyz))
            elif signum == 4:
                f.write("  {:<2} {:10.4f}{:12.4f}{:12.4f}\n".format(element, *xyz))
        f.write('end\n')
        f.write('dft\n')
        f.write('  xc b3lyp\n')
        f.write('  mult {}\n'.format(SpinNum))
        f.write('  iterations 200\n')
        f.write('end\n')
        f.write('task dft optimize\n')

def write_nwchem_fcf(nwcfcf2, smchg, SpinNum):
    """
    Writes an input file for a frequency calculation in NWChem.

    Parameters:
        filename (str): The name of the input file to write.
        smchg (int): The charge of the system.
        SpinNum (int): The multiplicity of the system.
        gatms (list): A list of `MoleculeAtom` objects representing the atoms in the system.
        signum (int): The number of decimal places to use for the atomic coordinates.
    """
    with open(nwcfcf2, 'w') as f:
        f.write('restart small_model\n')
        f.write('title "Frequency Calculation"\n')
        f.write('charge {}\n'.format(smchg))
        f.write('basis\n')
        f.write('  * library 6-31g*\n')
        f.write('end\n')
        f.write('dft\n')
        f.write('  xc b3lyp\n')
        f.write('  mult {}\n'.format(SpinNum))
        f.write('  iterations 200\n')
        f.write('end\n')
        f.write('task dft frequencies\n')


def write_nwchem_mkf(nwcmkf, lgchg, SpinNum, gatms, largeopt, signum=3):
    """
    Writes an input file for an ESP charge calculation in NWChem.

    Parameters:
        filename (str): The name of the input file to write.
        lgchg (int): The charge of the system.
        SpinNum (int): The multiplicity of the system.
        gatms (list): A list of `MoleculeAtom` objects representing the atoms in the system.
        signum (int): The number of decimal places to use for the atomic coordinates.
    """
    with open(nwcmkf, 'w') as f:
        f.write('start large_model\n')
        f.write('title "ESP Charge Calculation"\n')
        f.write('charge {}\n'.format(lgchg))
        f.write('basis\n')
        f.write('  * library 6-31g*\n')
        f.write('end\n')
        f.write('geometry\n')
        for gatm in gatms:
            element = gatm.element
            #nuchg = AtomicNum[element]
            #nuchg = round(nuchg, signum)
            xyz = gatm.crdx, gatm.crdy, gatm.crdz
            if signum == 3:
                f.write("  {:<2} {:10.3f}{:10.3f}{:10.3f}\n".format(element, *xyz))
            elif signum == 4:
                f.write("  {:<2} {:10.4f}{:12.4f}{:12.4f}\n".format(element, *xyz))
        f.write('end\n')
        f.write('dft\n')
        f.write('  xc b3lyp\n')
        f.write('  mult {}\n'.format(SpinNum))
        f.write('  iterations 200\n')
        f.write('end\n')
        if largeopt in [1, 2]:
            f.write('task dft optimize\n')
        # consider adding custom radius like in gauio
        f.write('esp\n')
        f.write('  recalculate\n')
        f.write('end\n')
        f.write('task ecp')

        