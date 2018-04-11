from unittest import TestCase
from atomium.structures import Model, Atom, Residue, Molecule
import biometal

class SolvationTests(TestCase):

    def setUp(self):
        self.model = Model()

        atom1 = Atom("Fe", 0, 0, 0, name="F1")
        mol = Molecule(atom1)
        self.model.add(mol)

        atom2 = Atom("C", 0.5, 0, 0, name="C")
        atom3 = Atom("S", 0, 0.5, 0, name="S1")
        atom4 = Atom("N", 0, 0, 0.5, name="N")
        atom5 = Atom("O", 0.5, 0, 0.5, name="O")
        res1 = Residue(atom2, atom3, atom4, atom5, name="VAL")
        self.model.add(res1)

        atom6 = Atom("C", 2.5, 0, 0, name="CA")
        atom7 = Atom("O", 0, 2.5, 0, name="O1", charge=-2)
        atom8 = Atom("N", 0, 0, 2.5, name="N1", charge=1)
        atom9 = Atom("O", 2, 0, 2, name="O")
        res2 = Residue(atom6, atom7, atom8, atom9, name="TRP")
        self.model.add(res2)

        atom10 = Atom("C", 4.5, 0, 0, name="CA")
        atom11 = Atom("O", 0, 4.5, 0, name="OE1")
        atom12 = Atom("O", 0, 0, 4.5, name="OE2")
        atom13 = Atom("O", 3, 0, 3, name="O2")
        res3 = Residue(atom10, atom11, atom12, atom13, name="GLU")
        self.model.add(res3)

        atom14 = Atom("C", 8.5, 0, 0, name="CP")
        atom15 = Atom("N", 0, 8.5, 0, name="ND1")
        atom16 = Atom("N", 0, 0, 8.5, name="NE2")
        atom17 = Atom("N", 6, 0, 6, name="N2")
        res4 = Residue(atom14, atom15, atom16, atom17, name="HIS")
        self.model.add(res4)


    def test_solvation_measures(self):
        self.assertEqual(biometal.solvation(self.model, 0, 0, 0, 0.1), 0)
        self.assertEqual(biometal.solvation(self.model, 0, 0, 0, 1), -1)
        self.assertEqual(biometal.solvation(self.model, 0, 0, 0, 2), -1)
        self.assertEqual(biometal.solvation(self.model, 0, 0, 0, 3), -71/9)
        self.assertEqual(biometal.solvation(self.model, 0, 0, 0, 4), -71/9)
        self.assertEqual(biometal.solvation(self.model, 0, 0, 0, 5), -108/13)
        self.assertEqual(biometal.solvation(self.model, 0, 0, 0, 6), -108/13)
        self.assertEqual(biometal.solvation(self.model, 0, 0, 0, 7), -108/13)
        self.assertEqual(biometal.solvation(self.model, 0, 0, 0, 8), -108/13)
        self.assertEqual(biometal.solvation(self.model, 0, 0, 0, 9), -146/17)
        self.assertEqual(
         biometal.solvation(self.model, 0, 0, 0, 9, het=False), -146/16
        )
        self.assertEqual(
         biometal.solvation(self.model, 0, 0, 0, 9, metal=False), -146/16
        )


    def test_partial_charges_measures(self):
        self.assertAlmostEqual(biometal.solvation(
         self.model, 0, 0, 0, 0.1, pc=True
        ), 0, delta=0.005)
        self.assertAlmostEqual(biometal.solvation(
         self.model, 0, 0, 0, 1, pc=True
        ), 0.797076/5, delta=0.005)
        self.assertAlmostEqual(biometal.solvation(
         self.model, 0, 0, 0, 2, pc=True
        ), 0.797076/5, delta=0.005)
        self.assertAlmostEqual(biometal.solvation(
         self.model, 0, 0, 0, 3, pc=True
        ), 6.10858/9, delta=0.005)
        self.assertAlmostEqual(biometal.solvation(
         self.model, 0, 0, 0, 4, pc=True
        ), 6.10858/9, delta=0.005)
        self.assertAlmostEqual(biometal.solvation(
         self.model, 0, 0, 0, 5, pc=True
        ), 7.165968/13, delta=0.005)
        self.assertAlmostEqual(biometal.solvation(
         self.model, 0, 0, 0, 6, pc=True
        ), 7.165968/13, delta=0.005)
        self.assertAlmostEqual(biometal.solvation(
         self.model, 0, 0, 0, 7, pc=True
        ), 7.165968/13, delta=0.005)
        self.assertAlmostEqual(biometal.solvation(
         self.model, 0, 0, 0, 8, pc=True
        ), 7.165968/13, delta=0.005)
        self.assertAlmostEqual(biometal.solvation(
         self.model, 0, 0, 0, 9, pc=True
        ), 7.688794/17, delta=0.005)



class ContrastTests(TestCase):

    def setUp(self):
        self.model = Model()

        atom1 = Atom("Fe", 0, 0, 0, name="F1")
        mol = Molecule(atom1)
        self.model.add(mol)

        atom2 = Atom("N", 0.5, 0, 0, name="N")
        atom3 = Atom("O", 0, 0.5, 0, name="O")
        res1 = Residue(atom2, atom3, name="TYR")
        self.model.add(res1)

        atom4 = Atom("C", 1, 0, 0, name="CA")
        atom5 = Atom("C", 0, 1, 0, name="CB")
        res2 = Residue(atom4, atom5, name="VAL")
        self.model.add(res2)


    def test_hydrophobic_contrast(self):
        # C = Σσr^2 - (n * σ r2)
        self.assertEqual(
         biometal.hydrophobic_contrast(self.model, 0, 0, 0, 0.1), 0
        )
        self.assertEqual(
         biometal.hydrophobic_contrast(self.model, 0, 0, 0, 0.5),
         ((0 * 0) + (0.25 * -9) + (0.25 * -9)) - (
          3 * ((0 + (-9) + (-9)) / 3) * ((0 + 0.25 + 0.25) / 3)
         )
        )
        self.assertEqual(
         biometal.hydrophobic_contrast(self.model, 0, 0, 0, 0.5, het=False),
         ((0.25 * -9) + (0.25 * -9)) - (2 * ((-9 + -9) / 2) * ((0.25 + 0.25) / 2))
        )
        self.assertEqual(
         biometal.hydrophobic_contrast(self.model, 0, 0, 0, 0.5, metal=False),
         ((0.25 * -9) + (0.25 * -9)) - (2 * ((-9 + -9) / 2) * ((0.25 + 0.25) / 2))
        )
        self.assertEqual(
         biometal.hydrophobic_contrast(self.model, 0, 0, 0, 1),
         ((0 * 0) + (0.25 * -9) + (0.25 * -9) + 18 + 18) - (
          5 * ((0 + (-9) + (-9) + 18 + 18) / 5) * ((0 + 0.25 + 0.25 + 1 + 1) / 5)
         )
        )
        self.assertEqual(
         biometal.hydrophobic_contrast(self.model, 0, 0, 0, 1, het=False),
         ((0.25 * -9) + (0.25 * -9) + 18 + 18) - (
          4 * (((-9) + (-9) + 18 + 18) / 4) * ((0.25 + 0.25 + 1 + 1) / 4)
         )
        )
        self.assertAlmostEqual(
         biometal.hydrophobic_contrast(self.model, 0, 0, 0, 1, pc=True),
         ((0 * 0) + (0.25 * -0.52 * -0.52) + (0.25 * -0.5 * -0.5) + (0.201 ** 2) + (0.033 ** 2)) - (
          5 * ((0 + ((-0.5) ** 2) + ((-0.52) ** 2) + (0.201 ** 2) + (0.033 ** 2)) / 5) * ((0 + 0.25 + 0.25 + 1 + 1) / 5)
         ), delta=0.0000005
        )
