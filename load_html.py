import csv
import os
import re
from pathlib import Path

from pyquery import PyQuery as pq

from logger import logger

htmls = sorted(list(Path("htmls").rglob("*.html")), key=lambda x: int(x.stem.split("_")[1]))

titles = []
urls = []
dois = []

for html_file in htmls:
    page = html_file.name.split(".")[0].split("_")[1]
    with open(html_file, "r", encoding='utf-8') as f:
        html = f.read()

    doc = pq(html)
    records = doc('app-record').items()

    count = 0
    for record in records:
        title = " ".join(record('app-summary-title a').text().split())
        url = record('app-summary-record-links a').attr("href")
        doi = re.findall("KeyAID=(.*)&DestApp=DOI", url) if url is not None else None
        doi = doi[0] if isinstance(doi, list) and len(doi) else None
        titles.append(title)
        urls.append(url)
        dois.append(doi)
        count += 1
    logger.info(f"Page {page}, record: {count}")

logger.info(f"Total titles: {len(titles)} urls: {len(urls)}")

data = list(zip(titles, urls, dois))

with open("datafile.csv", "w", encoding="utf-8", newline="") as f:
    head = ["title", 'url', 'doi']
    writer = csv.writer(f)
    writer.writerow(head)
    writer.writerows(data)

os.system("bash paper.sh datafile.csv> datafile.1.csv")
