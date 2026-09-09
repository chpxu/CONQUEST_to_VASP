
Transforms and supercells
=========================

C2a can handle transforming unit cells like in VESTA, with the important feature of being able to preserve species numbers. In CONQUEST, elements with different spins are given their own species, e.g. species 1 could be Mn (up) and species 2 could be Mn (down). Thus, creating and transforming cells with spin patterns is now trivial.

In VESTA, you can give atoms vectors. Transforming a unit cell preserves symmetry and thus the arrows (by choosing "equivalent sites"). However, once the symmetry is removed, the vectors are not maintained (understandably so). Since VESTA has no knowledge of "species", it has no notion that a set of atoms are meant to have particular spins, and that spins are mapped to species. For small cells this isn't an issue, you could set atoms to different elements temporarily and manually edit the coordinates file afterwards, but for large cells this is unfeasible and a time waste, which is where this module aims to fill the gap as CONQUEST is designed for large-scale simulations.

.. automodule:: conquest2a.cell.supercell
  :members:

.. automodule:: conquest2a.cell.transform
  :members: