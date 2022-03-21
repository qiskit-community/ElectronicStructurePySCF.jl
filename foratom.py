import electronic_structure as es
from electronic_structure import Atom, Geometry,  MolecularSpec
import numpy as np


def get_spec():
    h2 = [["H", (0.1, 0.2, 0.3)],  ("H", [1.1, 0.2, 3.3])]
    ms = MolecularSpec(h2, basis="sto3g", ensure=True)
    return ms


def trysome():
#    h2 = Geometry([Atom("H", (0.1, 0.2, 0.3)),  Atom("Ne", (1.1, 0.2, 3.3))])
#    h2 = [("H", (0.1, 0.2, 0.3)),  ("Ne", (1.1, 0.2, 3.3))]
#    h2 = [["H", (0.1, 0.2, 0.3)],  ("Ne", [1.1, 0.2, 3.3])]
    h2 = [["H", np.array([0.1, 0.2, 0.3])],  ("Ne", [1.1, 0.2, 3.3])]
    res = ensure_geom(h2)
    return res
