import csv
from pathlib import Path

from pyquery import PyQuery as pq

from logger import logger

htmls = sorted(list(Path("htmls").rglob("*.html")), key=lambda x: int(x.stem.split("_")[1]))

titles = []
urls = []

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
        titles.append(title)
        urls.append(url)
        count += 1
    logger.info(f"Page {page}, record: {count}")

logger.info(f"Total titles: {len(titles)} urls: {len(urls)}")

data = list(zip(titles, urls))

with open("datafile.csv", "w", encoding="utf-8") as f:
    head = ["title", 'url']
    writer = csv.writer(f)
    writer.writerow(head)
    writer.writerows(data)
