import sys, os, re
from openff.toolkit.topology import Molecule
from openmm.app import PDBFile, Modeller, ForceField
from openmmforcefields.generators import SystemGenerator
from openff.units.openmm import to_openmm
import openff.units as openff_units
from openmm import unit
from openff.interchange import Interchange

os.environ["INTERCHANGE_EXPERIMENTAL"] = "1"


protein_pdb = PDBFile('protein.pdb')
ligand_molecule = Molecule.from_file('ligand.sdf')
#ligand_molecule.name = "MOL"
#for residue in ligand_molecule.to_topology().to_openmm().residues():
#    residue.name = 'MOL'

modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
modeller.add(
    ligand_molecule.to_topology().to_openmm(),
    to_openmm(ligand_molecule.conformers[0])
)

ff_for_templates = ForceField(
    'amber/ff14SB.xml',
    'amber/tip3p_standard.xml',
    "amber/tip3p_HFE_multivalent.xml"
)

from openmmforcefields.generators import EspalomaTemplateGenerator
espaloma = EspalomaTemplateGenerator(molecules = ligand_molecule)
ff_for_templates.registerTemplateGenerator(espaloma.generator)

"""
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
smirnoff_generator = SMIRNOFFTemplateGenerator(
    molecules = ligand_molecule,
    force_field = 'openff-2.2.0.offxml'
)
ff_for_templates.registerTemplateGenerator(smirnoff_generator.generator)
"""

modeller.addSolvent(
    ff_for_templates, boxShape = "octahedron",
    model = 'tip3p', padding = 1.0 * unit.nanometers,
    positiveIon = 'Na+', negativeIon = 'Cl-', ionicStrength = 0.15 * unit.molar
)


"""
# トポロジーと位置情報をそれぞれ取得
top = modeller.getTopology()
pos = modeller.getPositions()
PDBFile.writeFile(top, pos, open("test.pdb", "w"))
"""
