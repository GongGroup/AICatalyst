import re
from pathlib import Path
from pyquery import PyQuery as pq

from logger import logger

literatures = Path("literature").rglob("*.html")

logger.info(f"Check download literature completeness")
for literature in literatures:
    with open(literature, "r", encoding="utf-8") as f:
        html = f.read()
    doc = pq(html)
    conclusion = re.findall("[Cc]onclusion|[Ss]ummary", doc.text())
    if not len(conclusion):
        logger.info(f"{literature}: Search keyword `conclusion` failed, please check")
