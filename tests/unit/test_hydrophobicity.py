from atomium.structures import Model
from unittest import TestCase
from unittest.mock import Mock, patch, MagicMock
from biometal.hydrophobicity import solvation, atom_solvation

class SolvationTests(TestCase):

    def setUp(self):
        self.model = Mock(Model)
        self.atoms = [Mock(), Mock(), Mock(), Mock(), Mock()]
        self.model.atoms.return_value = self.atoms
        self.patch1 = patch("biometal.hydrophobicity.Atom")
        self.patch2 = patch("biometal.hydrophobicity.atom_solvation")
        self.mock_atom = self.patch1.start()
        self.mock_atsolv = self.patch2.start()
        self.dummy = Mock()
        self.dummy.nearby.return_value = self.atoms[:3]
        self.mock_atom.return_value = self.dummy


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
        self.mock_atom.assert_called_with("X", 2, 4, 5)
        self.model.add_atom.assert_called_with(self.dummy)
        self.dummy.nearby.assert_called_with(12)
        self.mock_atsolv.assert_any_call(self.atoms[0])
        self.mock_atsolv.assert_any_call(self.atoms[1])
        self.mock_atsolv.assert_any_call(self.atoms[2])
        self.model.remove_atom.assert_called_with(self.dummy)
        self.assertEqual(solv, 6)


    def test_dummy_atom_is_always_removed(self):
        self.dummy.nearby.side_effect = Exception
        with self.assertRaises(Exception):
            solvation(self.model, 2, 4, 5, 12)
        self.model.remove_atom.assert_called_with(self.dummy)



class AtomSolvationTests(TestCase):

    def setUp(self):
        self.atom = Mock()


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
