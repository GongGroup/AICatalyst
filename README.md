# AICatalyst Manual

## Table of Contents

- [About AICatalyst](#about-aicatalyst)
- [Code Structure](#code-structure)
- [Requirements](#requirements)

## About AICatalyst

This package mainly services for obtaining literatures from [Web of Science](https://www.webofscience.com) (WOS),
construct the database, as well as use the machine learning to design the catalyst.

## Code Structure

- [common](AICatalysis/common): common API

- [crawler](AICatalysis/crawler): carry the crawler task from WOS or publisher

- [database](AICatalysis/database): related code about database construction

- [tutorial](tutorial): a brief tutorial about the software and knowledge prepare to quickly start this project

## Requirements

- lxml
- rdkit
- numpy
- execjs
- pyquery
- selenium
- chromedriver
- requests_html
- Beautiful Soup
- undetected_chromedriver
