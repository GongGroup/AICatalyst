from pathlib import Path

from common.fio import JsonIO

# Directory constant here
ChemDir = Path("../chemical")
DrawDir = Path("../draw")

# file constant
FChemical = ChemDir / "chemical.json"
FReaxys = ChemDir / "reaxys.json"
FReactions = ChemDir / "reactions.json"
FReaxysYield = ChemDir / "opsin_reaxys.json"
FOpsinRecord = ChemDir / "opsin_record.json"
FRXconfig = ChemDir / "RX_config.json"
FRXDconfig = ChemDir / "RXD_config.json"
FReaxysXML = ChemDir / "reaxys_xml.xml"
FSucReaxys = ChemDir / "suc_reaxys.json"

# Variable constant
ChemInfo = {item['name']: {key: value for key, value in item.items() if key != 'name'}
            for item in JsonIO.read(FOpsinRecord)}
