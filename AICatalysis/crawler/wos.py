import csv
import hashlib
import json
import logging
import random
import re
import time
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from pyquery import PyQuery
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as ec

from AICatalysis.common.driver import ChromeDriver
from AICatalysis.common.logger import init_root_logger

WOSRoot = "https://www.webofscience.com/wos/alldb/basic-search"
DataDir = Path("../../data")

init_root_logger()
logger = logging.getLogger(__name__)


class WOSCrawler(object):
    """
    According to the setting in config.json, crawler the WOS database and write the results to output/ directory,
    one html represents one page in WOS.
    """

    def __init__(self, config_file, output):
        """
        Initialize the WOSCrawler

        Args:
            config_file: configuration path in json format
            output: represent the output path
        """
        self.start_page, self.end_page = None, None
        self.keyword = None
        self.configure(config_file)
        self.output = output

    def configure(self, config_file):
        with open(f"{config_file}", "r") as f:
            config = json.load(f)

        # load config setting
        self.start_page = config['start_page']
        self.end_page = config['end_page']
        self.keyword = config['keyword']

    def get_htmls(self):
        # start driver

        chrome = ChromeDriver(timeout=60, options=False)
        driver, wait = chrome.driver, chrome.wait
        driver.get(WOSRoot)

        logger.info("prepare to accept cookie")
        cookie_accept = wait.until(ec.visibility_of_element_located((By.CSS_SELECTOR, '#onetrust-accept-btn-handler')))
        driver.execute_script("arguments[0].click();", cookie_accept)

        logger.info("prepare to input search")
        search_input = wait.until(ec.visibility_of_element_located((By.XPATH, '//*[@id="mat-input-0"]')))
        search_input.click()
        search_input.send_keys(self.keyword)

        logger.info("prepare to submit")
        submit_button = wait.until(
            ec.element_to_be_clickable((By.XPATH, '//*[@id="snSearchType"]/div[3]/button[2]/span[1]')))
        submit_button.click()
        time.sleep(random.random() * 9)

        # obtain total pages
        cur_url = driver.current_url
        time.sleep(random.random())
        doc = PyQuery(driver.page_source)
        tot_pages = int(doc("span[class=end-page\ ng-star-inserted]").text().split()[0])
        logger.info(f"Total Pages: {tot_pages}")

        if self.start_page > tot_pages:
            raise ValueError(f"start_page = {self.start_page} over total_pages({tot_pages})")

        if self.end_page > tot_pages:
            logger.warning(f"end_page = {self.end_page} over total_pages({tot_pages}), reset end_page = {tot_pages}")
            self.end_page = tot_pages

        cur_page = self.start_page

        # start downloading
        while cur_page <= self.end_page:

            # set random sleep
            sleep_t1 = random.random() + 1
            sleep_t3 = sleep_t1 * 3
            sleep_t5 = sleep_t1 * 5
            sleep_t9 = sleep_t1 * 9

            logger.info(f"Page {cur_page} downloading...")
            cur_url = "/".join(cur_url.split("/")[:-1]) + f"/{cur_page}"
            driver.get(cur_url)
            time.sleep(sleep_t3)

            # obtain the scrollHeight by js
            scroll_height = driver.execute_script("var q=document.body.scrollHeight; return(q)")

            # generate random scroll distance, total scroll_t times, sum is scrollHeight
            scroll_t = random.randint(20, 40)
            loop_height = np.random.dirichlet(np.ones(scroll_t)) * scroll_height
            for i, h in zip(range(scroll_t), loop_height):
                driver.execute_script(f"window.scrollBy(0,{h})")
                time.sleep(sleep_t1)
            driver.execute_script(f"window.scrollTo(0,{scroll_height})")

            # write to *.html files
            with open(f"{self.output}/page_{cur_page}.html", "w", encoding='utf-8') as f:
                f.write(driver.page_source)
            cur_page += 1

            if cur_page % 5 == 0:
                time.sleep(sleep_t5)

            if cur_page % 20 == 0:
                time.sleep(sleep_t9)

        logger.info("All Pages download successful")

        driver.close()


class WOSParser(object):
    """
    Load database from htmls/ directory, parse title, url and doi field, output to a csv file
    """

    def __init__(self, htmls_dir: str, datafile: Union[str, Path], with_url: bool = True, md5_flag: bool = True,
                 jcr: bool = True):
        """
        Initialize the WOSCrawler

        Args:
            htmls_dir (str): htmls directory store wos database
            datafile (str): output path in csv format
            with_url (bool): whether output record without url
        """
        self.htmls_dir = sorted(list(Path(htmls_dir).rglob("*.html")), key=lambda x: int(x.stem.split("_")[1]))
        self.datafile = datafile
        self.with_url = with_url
        self.md5_flag = md5_flag
        if jcr:
            self.jcr = pd.read_excel(DataDir / "2021JCR.xlsx")

    def __call__(self):
        self.parse()

    def parse(self):

        titles = []
        urls = []
        dois = []
        years = []
        journals = []
        if_ = []

        for html_file in self.htmls_dir:
            page = html_file.name.split(".")[0].split("_")[1]
            with open(html_file, "r", encoding='utf-8') as f:
                html = f.read()

            doc = PyQuery(html)
            records = doc('app-record').items()

            tot_count, effective_count = 0, 0
            for record in records:
                tot_count += 1
                title = " ".join(record('app-summary-title a').text().split())
                url = record('app-summary-record-links a').attr("href")
                if self.with_url and url is None:
                    continue
                doi = re.findall("KeyAID=(.*)&DestApp=DOI", url) if url is not None else None
                doi = doi[0] if isinstance(doi, list) and len(doi) else None

                # parse publish year
                year = record('span[name]').text()
                if not year:
                    year = record('span.value').text()  # e.g., Jun 2022 (在线发表) |
                year_re = re.compile("\d{4}")
                fin_year = ""
                for item in re.split(" |-|\)", year):
                    if re.match(year_re, item):
                        fin_year = item
                        break
                year = int(fin_year)

                # parse journal && IF
                journal = record('app-jcr-overlay').text()
                try:
                    if_2022 = self.jcr[self.jcr['journal_name'].str.fullmatch(journal, case=False)]['if_2022'].values[0]
                except IndexError:
                    if_2022 = 0.

                titles.append(title)
                urls.append(url)
                dois.append(doi)
                years.append(year)
                journals.append(journal)
                if_.append(if_2022)

                effective_count += 1

            logger.info(f"Page {page}, record_tot: {tot_count}, record_effective: {effective_count}")

        logger.info(f"Total titles: {len(titles)} urls: {len(urls)}")
        exit()
        data = list(zip(titles, urls, dois))

        if self.md5_flag:
            md5_dict = {hashlib.md5(url.encode(encoding='utf-8')).hexdigest(): url for title, url, doi in data}
            with open("../../data/md5_name.json", "w", encoding="utf-8") as f:
                json.dump(md5_dict, f)
            logger.info("url~md5 mapping has been wrote to md5_name.json")

        with open(f"{self.datafile}", "w", encoding="utf-8", newline="") as f:
            head = ["title", 'url', 'doi']
            writer = csv.writer(f)
            writer.writerow(head)
            writer.writerows(data)


if __name__ == '__main__':
    # crawler = WOSCrawler(config_file="config.json", output="htmls")
    # crawler.get_htmls()

    parser = WOSParser(htmls_dir="../../htmls", datafile=DataDir / "datafile.csv")
    parser()
