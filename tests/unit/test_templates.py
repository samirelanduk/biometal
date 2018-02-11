from atomium.structures.chains import Site
from atomium.structures.molecules import AtomicStructure
from unittest import TestCase
from unittest.mock import Mock, MagicMock, patch
from biometal.templates import create_site_template

class TemplateCreationTests(TestCase):

    def setUp(self):
        self.site = Mock(Site)
        self.residues = [Mock(), Mock(), Mock()]
        self.atoms = ["A{}".format(i) for i in range(1, 21)]
        self.residues[0].atom.side_effect = self.atoms[:2]
        self.residues[1].atom.side_effect = self.atoms[2:4]
        self.residues[2].atom.side_effect = self.atoms[4:6]
        self.site.residues.return_value = self.residues
        self.patch1 = patch("biometal.templates.AtomicStructure")
        self.mock_struct = self.patch1.start()


    def tearDown(self):
        self.patch1.stop()


    def test_template_creation_needs_site(self):
        with self.assertRaises(TypeError):
            create_site_template("site")


    def test_can_create_template(self):
        template = create_site_template(self.site)
        self.mock_struct.assert_called_with(*self.atoms[:6])
        for res in self.residues:
            res.atom.assert_any_call(name="CA")
            res.atom.assert_any_call(name="CB")
