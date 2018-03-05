from atomium.structures import Model
from unittest import TestCase
from unittest.mock import Mock, patch, MagicMock
from biometal.hydrophobicity import *

class SolvationTests(TestCase):

    def setUp(self):
        self.model = Mock(Model)
        self.atoms = [Mock(), Mock(), Mock(), Mock(), Mock()]
        self.model.atoms_in_sphere.return_value = self.atoms[:3]
        self.patch1 = patch("biometal.hydrophobicity.atom_solvation")
        self.patch2 = patch("biometal.hydrophobicity.atom_partial_charge")
        self.mock_atsolv = self.patch1.start()
        self.mock_atcharge = self.patch2.start()


    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()


    def test_solvation_needs_model(self):
        with self.assertRaises(TypeError):
            solvation("structure", 0, 0, 0, 10)


    def test_coordinates_must_be_numbers(self):
        with self.assertRaises(TypeError):
            solvation(self.model, "0", 0, 0, 10)
        with self.assertRaises(TypeError):
            solvation(self.model, 0, "0", 0, 10)
        with self.assertRaises(TypeError):
            solvation(self.model, 0, 0, "0", 10)


    def test_cutoff_must_be_positive_number(self):
        with self.assertRaises(TypeError):
            solvation(self.model, 0, 0, 0, "10")
        with self.assertRaises(ValueError):
            solvation(self.model, 0, 0, 0, -10)


    def test_can_get_solvation(self):
        self.mock_atsolv.side_effect = [11, -9, 4]
        solv = solvation(self.model, 2, 4, 5, 12)
        self.model.atoms_in_sphere.assert_called_with(
         2, 4, 5, 12, het=True, metal=True
        )
        self.mock_atsolv.assert_any_call(self.atoms[0])
        self.mock_atsolv.assert_any_call(self.atoms[1])
        self.mock_atsolv.assert_any_call(self.atoms[2])
        self.assertFalse(self.mock_atcharge.called)
        self.assertEqual(solv, 2)


    def test_can_get_partial_charges(self):
        self.mock_atcharge.side_effect = [11, -8, 5]
        solv = solvation(self.model, 2, 4, 5, 12, pc=True)
        self.model.atoms_in_sphere.assert_called_with(
         2, 4, 5, 12, het=True, metal=True
        )
        self.mock_atcharge.assert_any_call(self.atoms[0])
        self.mock_atcharge.assert_any_call(self.atoms[1])
        self.mock_atcharge.assert_any_call(self.atoms[2])
        self.assertFalse(self.mock_atsolv.called)
        self.assertEqual(solv, 70)


    def test_can_handle_zero_atoms(self):
        self.model.atoms_in_sphere.return_value = []
        self.assertEqual(solvation(self.model, 2, 4, 5, 12), 0)


    def test_can_filter_out_heteroatoms(self):
        solvation(self.model, 2, 4, 5, 12, het=False)
        self.model.atoms_in_sphere.assert_called_with(
         2, 4, 5, 12, het=False, metal=True
        )


    def test_can_filter_out_metal(self):
        solvation(self.model, 2, 4, 5, 12, metal=False)
        self.model.atoms_in_sphere.assert_called_with(
         2, 4, 5, 12, het=True, metal=False
        )



class AtomSolvationTests(TestCase):

    def setUp(self):
        self.atom, self.residue = Mock(), Mock()
        self.atom.charge.return_value = 0
        self.atom.name.return_value = None
        self.atom.residue.return_value = self.residue
        self.residue.name.return_value = "MET"


    def test_default_solvation_is_zero(self):
        self.assertEqual(atom_solvation(self.atom), 0)


    def test_carbon_solvation(self):
        self.atom.element.return_value = "C"
        self.assertEqual(atom_solvation(self.atom), 18)


    def test_oxygen_solvation(self):
        self.atom.element.return_value = "O"
        self.assertEqual(atom_solvation(self.atom), -9)


    def test_nitrogen_solvation(self):
        self.atom.element.return_value = "N"
        self.assertEqual(atom_solvation(self.atom), -9)


    def test_charged_oxygen_solvation(self):
        self.atom.element.return_value = "O"
        self.atom.charge.return_value = -2
        self.assertEqual(atom_solvation(self.atom), -37)


    def test_charged_nitrogen_solvation(self):
        self.atom.element.return_value = "N"
        self.atom.charge.return_value = 1
        self.assertEqual(atom_solvation(self.atom), -38)


    def test_sulphur_solvation(self):
        self.atom.element.return_value = "S"
        self.assertEqual(atom_solvation(self.atom), -5)


    def test_glutamate_oxygens(self):
        self.residue.name.return_value = "GLU"
        self.atom.element.return_value = "O"
        self.assertEqual(atom_solvation(self.atom), -9)
        self.atom.name.return_value = "OE1"
        self.assertEqual(atom_solvation(self.atom), -23)
        self.atom.name.return_value = "OE2"
        self.assertEqual(atom_solvation(self.atom), -23)
        self.atom.name.return_value = "OE3"
        self.assertEqual(atom_solvation(self.atom), -9)
        self.atom.residue.return_value = None
        self.assertEqual(atom_solvation(self.atom), -9)


    def test_aspartate_oxygens(self):
        self.residue.name.return_value = "ASP"
        self.atom.element.return_value = "O"
        self.assertEqual(atom_solvation(self.atom), -9)
        self.atom.name.return_value = "OD1"
        self.assertEqual(atom_solvation(self.atom), -23)
        self.atom.name.return_value = "OD2"
        self.assertEqual(atom_solvation(self.atom), -23)
        self.atom.name.return_value = "OD3"
        self.assertEqual(atom_solvation(self.atom), -9)
        self.atom.residue.return_value = None
        self.assertEqual(atom_solvation(self.atom), -9)


    def test_histidine_nitrogens(self):
        self.residue.name.return_value = "HIS"
        self.atom.element.return_value = "N"
        self.assertEqual(atom_solvation(self.atom), -9)
        self.atom.name.return_value = "ND1"
        self.assertEqual(atom_solvation(self.atom), -23.5)
        self.atom.name.return_value = "NE2"
        self.assertEqual(atom_solvation(self.atom), -23.5)
        self.atom.name.return_value = "NF3"
        self.assertEqual(atom_solvation(self.atom), -9)
        self.atom.residue.return_value = None
        self.assertEqual(atom_solvation(self.atom), -9)


    def test_arginine_nitrogens(self):
        self.residue.name.return_value = "ARG"
        self.atom.element.return_value = "N"
        self.assertEqual(atom_solvation(self.atom), -9)
        self.atom.name.return_value = "NH1"
        self.assertEqual(atom_solvation(self.atom), -23.5)
        self.atom.name.return_value = "NH2"
        self.assertEqual(atom_solvation(self.atom), -23.5)
        self.atom.name.return_value = "NH3"
        self.assertEqual(atom_solvation(self.atom), -9)
        self.atom.residue.return_value = None
        self.assertEqual(atom_solvation(self.atom), -9)



class AtomicPartialChargesTests(TestCase):

    def setUp(self):
        self.atom, self.residue = Mock(), Mock()
        self.atom.charge.return_value = 0
        self.atom.name.return_value = None
        self.atom.residue.return_value = self.residue
        self.patch1 = patch("biometal.hydrophobicity.partial_charges")
        self.mock_charges = self.patch1.start()
        d = {"AAA": {"N": 1, "M": 2}, "BBB": {"X": 3, "Y": 4}}
        self.mock_charges.__getitem__.side_effect = d.__getitem__
        self.mock_charges.__contains__.side_effect = d.__contains__


    def tearDown(self):
        self.patch1.stop()


    def test_default_solvation_is_zero(self):
        self.assertEqual(atom_partial_charge(self.atom), 0)


    def test_can_return_actual_charge(self):
        self.atom.charge.return_value = -2
        self.assertEqual(atom_partial_charge(self.atom), -2)
        self.atom.charge.return_value = 1.5
        self.assertEqual(atom_partial_charge(self.atom), 1.5)


    def test_can_return_partial_charge(self):
        self.residue.name.return_value = "AAA"
        self.atom.name.return_value = "N"
        self.assertEqual(atom_partial_charge(self.atom), 1)
        self.atom.name.return_value = "M"
        self.assertEqual(atom_partial_charge(self.atom), 2)
        self.residue.name.return_value = "BBB"
        self.atom.name.return_value = "X"
        self.assertEqual(atom_partial_charge(self.atom), 3)
        self.atom.name.return_value = "Y"
        self.assertEqual(atom_partial_charge(self.atom), 4)


    def test_can_handle_unknown_atoms(self):
        self.residue.name.return_value = "AAA"
        self.atom.name.return_value = "N"
        self.assertEqual(atom_partial_charge(self.atom), 1)
        self.atom.name.return_value = "P"
        self.assertEqual(atom_partial_charge(self.atom), 0)
        self.residue.name.return_value = "CCC"
        self.assertEqual(atom_partial_charge(self.atom), 0)
        self.atom.residue.return_value = None
        self.assertEqual(atom_partial_charge(self.atom), 0)



class HydrophobicContrastTests(TestCase):

    def setUp(self):
        self.model = Mock(Model)
        self.atoms = [Mock(), Mock(), Mock(), Mock(), Mock()]
        self.atoms[0].distance_to.return_value = 7
        self.atoms[1].distance_to.return_value = 5
        self.atoms[2].distance_to.return_value = 10
        self.model.atoms_in_sphere.return_value = self.atoms[:3]
        self.patch1 = patch("biometal.hydrophobicity.solvation")
        self.patch2 = patch("biometal.hydrophobicity.atom_solvation")
        self.mock_solv = self.patch1.start()
        self.mock_atsolv = self.patch2.start()
        self.mock_solv.return_value = 8


    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()


    def test_contrast_needs_model(self):
        with self.assertRaises(TypeError):
            hydrophobic_contrast("structure", 0, 0, 0, 10)


    def test_coordinates_must_be_numbers(self):
        with self.assertRaises(TypeError):
            hydrophobic_contrast(self.model, "0", 0, 0, 10)
        with self.assertRaises(TypeError):
            hydrophobic_contrast(self.model, 0, "0", 0, 10)
        with self.assertRaises(TypeError):
            hydrophobic_contrast(self.model, 0, 0, "0", 10)


    def test_cutoff_must_be_positive_number(self):
        with self.assertRaises(TypeError):
            hydrophobic_contrast(self.model, 0, 0, 0, "10")
        with self.assertRaises(ValueError):
            hydrophobic_contrast(self.model, 0, 0, 0, -10)


    def test_can_get_hydrophobic_contrast(self):
        self.mock_atsolv.side_effect = [11, -9, 4]
        contrast = hydrophobic_contrast(self.model, 4, 8, 15, 10)
        self.model.atoms_in_sphere.assert_called_with(
         4, 8, 15, 10, het=True, metal=True
        )
        self.mock_solv.assert_called_with(
         self.model, 4, 8, 15, 10, het=True, metal=True
        )
        for atom in self.atoms[:3]:
            atom.distance_to.assert_any_call((4, 8, 15))
            self.mock_atsolv.assert_any_call(atom)
        self.assertEqual(
         contrast, ((11 * 49) + (-9 * 25) + (4 * 100)) - (3 * 8 * 58)
        )


    def test_can_handle_zero_atoms(self):
        self.model.atoms_in_sphere.return_value = []
        self.assertEqual(hydrophobic_contrast(self.model, 2, 4, 5, 12), 0)


    def test_can_filter_out_heteroatoms(self):
        hydrophobic_contrast(self.model, 2, 4, 5, 12, het=False)
        self.model.atoms_in_sphere.assert_called_with(
         2, 4, 5, 12, het=False, metal=True
        )
        self.mock_solv.assert_called_with(
         self.model, 2, 4, 5, 12, het=False, metal=True
        )


    def test_can_filter_out_metal(self):
        hydrophobic_contrast(self.model, 2, 4, 5, 12, metal=False)
        self.model.atoms_in_sphere.assert_called_with(
         2, 4, 5, 12, het=True, metal=False
        )
        self.mock_solv.assert_called_with(
         self.model, 2, 4, 5, 12, het=True, metal=False
        )
