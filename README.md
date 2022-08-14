# WOSCrawler Manual

## Table of Contents

- [About WOSCrawler](#about-woscrawler)
- [Code Structure](#code-structure)
- [Requirements](#requirements)

## About WOSCrawler

This package mainly services for obtaining literatures from [Web of Science](https://www.webofscience.com) (WOS).

## Code Structure

- [driver](driver.py): customized selenium-driver

- [wos](wos.py): crawler and parse `*.html` from WOS

- [publisher](publisher.py): according the url to crawler and check paper in `html` format

- [parse_reaxys](parse_reaxys.py): in charge to parse [reaxys](https://www.reaxys.com/) database

- [ichem](ichem.py): obtain structure in `sdf` format according to chemical name

## Requirements

- numpy
- pyquery
- selenium
- chromedriver
- Beautiful Soup
- undetected_chromedriver