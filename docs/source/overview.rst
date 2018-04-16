Overview
--------

biometal analyses metalloproteins.

Hydrophobic Contrast
~~~~~~~~~~~~~~~~~~~~

biometal can implement the hydrophobic contrast function as implemented in
(1990) Where metal ions bind in proteins:

  >>> import atomium
  >>> import biometal
  >>> model = atomium.fetch('1ton').model
  >>> zinc = model.atom(element='ZN')
  >>> biometal.hydrophobic_contrast(model, *zinc.location, 4)
  435.9691187500001

This uses atomic solvation parameters to determine hydrophobicity. To use
the square of partial charge instead, use the ``pc=True`` argument.
