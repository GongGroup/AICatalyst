import json
import re
from pathlib import Path

from lxml.etree import ParserError
from pyquery import PyQuery as pq

from logger import logger

# load excluded
with open("exclude", "r", encoding="utf-8") as f:
    exclude = [item.split(",")[0] for item in f.readlines()]

# load checked
with open("checked", "r", encoding="utf-8") as f:
    checked = [item.strip() for item in f.readlines()]

# load md5_dict
with open("md5_name.json", "r", encoding="utf-8") as f:
    md5_dict = json.load(f)

# search which html failed
literatures = Path("literature").rglob("*.html")
logger.info(f"Check download literature completeness")
for literature in literatures:
    if literature.stem in checked:
        continue
    conclusion = []
    with open(literature, "r", encoding="utf-8") as f:
        html = f.read()
    try:
        doc = pq(html)
        conclusion = re.findall("[Cc]onclusion|[Ss]ummary|[Oo]utlook|[Dd]iscussion|[Ee]xperiment|总结|结论", doc.text())
    except ParserError:
        pass

    if not len(conclusion) and literature.stem not in exclude:
        logger.info(f"{literature}: Search keyword `conclusion` failed, prepare check")
        logger.info(md5_dict[literature.stem.strip()])
    elif len(conclusion):
        with open("checked", "a+", encoding="utf-8") as f:
            f.write(f"{literature.stem} \n")
        logger.info(f"{literature} store in database")
