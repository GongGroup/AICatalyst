## WOSCrawlel

This package mainly services for obtaining literatures from [`Web of Science`](https://www.webofscience.com)

- get_html: obtain the *.html by **selenium**

- load_html: parse *.html by **pyquery** and generate the datafile.csv including `title` and `url` fields

- download: download literature from url

- download_check: check effectivity of download literature

## Requirements

- numpy
- pyquery
- selenium
- chromedriver