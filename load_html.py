import csv
from pathlib import Path

from pyquery import PyQuery as pq

from logger import logger

htmls = Path("htmls").rglob("*.html")

titles = []
urls = []

for html_file in htmls:
    page = html_file.name.split(".")[0].split("_")[1]
    with open(html_file, "r", encoding='utf-8') as f:
        html = f.read()

    doc = pq(html)
    records = doc('app-record').items()

    for record in records:
        title = " ".join(record('app-summary-title a').text().split())
        url = record('app-summary-record-links a').attr("href")
        titles.append(title)
        urls.append(url)

logger.info(f"Total titles: {len(titles)} urls: {len(urls)}")

data = list(zip(titles, urls))

with open("datafile.csv", "w", encoding="utf-8") as f:
    head = ["title", 'url']
    writer = csv.writer(f)
    writer.writerow(head)
    writer.writerows(data)
