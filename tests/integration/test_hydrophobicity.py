from unittest import TestCase
from atomium.structures import Model, Atom, Residue, Molecule
import biometal

class Tests(TestCase):

    def setUp(self):
        self.model = Model()

        atom1 = Atom("Fe", 0, 0, 0, name="F1")
        mol = Molecule(atom1)
        self.model.add_molecule(mol)

        atom2 = Atom("C", 0.5, 0, 0, name="CA")
        atom3 = Atom("C", 0, 0.5, 0, name="CB")
        atom4 = Atom("N", 0, 0, 0.5, name="N1")
        atom5 = Atom("C", 0.5, 0, 0.5, name="CC")
        res1 = Residue(atom2, atom3, atom4, atom5, name="VAL")
        self.model.add_residue(res1)


    def test_solvation_measures(self):
        self.assertEqual(biometal.solvation(self.model, 0, 0, 0, 0.1), 0)
        self.assertEqual(biometal.solvation(self.model, 0, 0, 0, 1), 45)
