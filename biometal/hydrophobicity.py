from atomium.structures import Model, Atom

def solvation(model, x, y, z, radius):
    """Determines the average solvation within a given sphere of an atomium
    model.

    :param Model model: The atomium model to examine.
    :param x: The x-coordinate of the centre of the sphere.
    :param y: The y-coordinate of the centre of the sphere.
    :param z: The z-coordinate of the centre of the sphere.
    :param radius: The radius of the sphere.
    :raises TypeError: if the model is not an atomium model object.
    :raises TypeError: if the coordinates are not numeric.
    :raises TypeError: if the radius is not numeric.
    :raises ValueError: if the radius is negative.
    :rtype: ``float``"""

    if not isinstance(model, Model):
        raise TypeError("{} is not a Model".format(model))
    if any(not isinstance(c, (int, float)) for c in (x, y, z)):
        raise TypeError("({}, {}, {}) not valid coordinate".format(x, y, z))
    if not isinstance(radius, (int, float)):
        raise TypeError("{} is not a valid radius".format(radius))
    if radius < 0:
        raise ValueError("{} is not a valid radius".format(radius))
    dummy = Atom("X", x, y, z)
    model.add_atom(dummy)
    try:
        sphere = dummy.nearby(radius)
        solvations = [atom_solvation(atom) for atom in sphere]
    finally:
        model.remove_atom(dummy)
    return sum(solvations)
    

def atom_solvation(atom):
    """Returns the atomic solvation parameter of an atomium atom.

    :param Atom atom: an atomium atom object.
    :rtype: ``float``"""

    if atom.element() == "C":
        return 18
    if atom.element() == "O" or atom.element() == "N":
        return -9
    return 0
