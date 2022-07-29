import json
import os
import re
from pathlib import Path

from pyquery import PyQuery as pq

from logger import logger

# load excluded
with open("exclude", "r", encoding="utf-8") as f:
    exclude = [item.split(",")[0] for item in f.readlines()]

# search which html failed
failed = []
literatures = Path("literature").rglob("*.html")
logger.info(f"Check download literature completeness")
for literature in literatures:
    with open(literature, "r", encoding="utf-8") as f:
        html = f.read()
    doc = pq(html)
    conclusion = re.findall("[Cc]onclusion|[Ss]ummary|[Oo]utlook|[Dd]iscussion|[Ee]xperiment|总结", doc.text())
    if not len(conclusion) and literature.stem not in exclude:
        logger.info(f"{literature}: Search keyword `conclusion` failed, prepare check")
        failed.append(literature.stem)

# load md5_dict
with open("md5_name.json", "r", encoding="utf-8") as f:
    md5_dict = json.load(f)

# locate which url failed
for name in failed:
    logger.info(md5_dict[name.strip()])
