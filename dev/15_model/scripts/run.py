from pathlib import Path

from libtbx.utils import null_out
from mmtbx import monomer_library
import mmtbx.model
import iotbx.pdb

if __name__ == '__main__':
    model = mmtbx.model.manager(
        model_input = iotbx.pdb.input(file_name=str(Path(Path.home(), "xray/data/pdbs/3ca7/3ca7.pdb"))),
        log=null_out()
    )
