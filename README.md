# WOSCrawler Manual

## About WOSCrawler

This package mainly services for obtaining literatures from [Web of Science](https://www.webofscience.com) (WOS).

## Code Structure

- [driver](driver.py): customized selenium-driver

- [wos](wos.py): crawler and parse `*.html` from WOS

- [publisher](publisher.py): according the url to crawler and check paper in `html` format

- [parse_reaxys](parse_reaxys.py): in charge to parse [reaxys](https://www.reaxys.com/) database

## Requirements

- numpy
- pyquery
- selenium
- chromedriver
- Beautiful Soup
- undetected_chromedriver