biometal
========

biometal is a suite of tools for dealing with metalloproteins.

Example
-------

  >>> import atomium
  >>> import biometal
  >>> model = atomium.fetch('1ton').model
  >>> zinc = model.atom(element='ZN')
  >>> biometal.hydrophobic_contrast(model, *zinc.location, 4)
  435.9691187500001



Table of Contents
-----------------

.. toctree ::

    installing
    api
    changelog
