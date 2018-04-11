from unittest import TestCase
import atomium
import biometal

class Tests(TestCase):

    def test_templates(self):
        site = atomium.fetch("1TON").model.molecule(name="ZN").site()

        template = biometal.create_site_template(site)
        self.assertEqual(len(template.atoms()), 6)
        self.assertEqual(len(template.atoms(name="CA")), 3)
        self.assertEqual(len(template.atoms(name="CB")), 3)
        template.save("test.pdb")
