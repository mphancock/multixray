import IMP
import IMP.atom

from pathlib import Path

m = IMP.Model()
hs = IMP.atom.read_multimodel_pdb(str(Path(Path.home(), "Documents/xray/dev/26_phenix_refine/data/test.pdb")), m, IMP.atom.NonWaterPDBSelector())
