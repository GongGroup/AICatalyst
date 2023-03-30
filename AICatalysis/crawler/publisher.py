import hashlib
import json
import logging
import random
import re
import time
from pathlib import Path
from typing import Union

import numpy as np
from lxml.etree import ParserError
from pyquery import PyQuery
from selenium.common import JavascriptException
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as ec

from AICatalysis.common.driver import ChromeDriver
from AICatalysis.common.logger import init_root_logger

DataDir = Path("../../data")

init_root_logger()
logger = logging.getLogger(__name__)


class PublisherCrawler(object):
    """
    According to url/doi in index_file, crawler corresponding publisher and stored them in output/ directory
    """

    def __init__(self, index_file: Union[str, Path], output):
        """
        Initialize the PublisherCrawler

        Args:
            index_file: stored the title, url and doi
            output: represent the output path
        """
        self.index_file = index_file
        self.output = output

    def get_htmls(self):
        with open(f"{self.index_file}", encoding='utf-8') as f:
            content = f.readlines()

        chrome = ChromeDriver()
        driver, wait = chrome.driver, chrome.wait
        scroll_height = 1000

        for index, line in enumerate(content[:5000]):
            url, doi = line.strip().split(",")[-2:]
            md5_name = hashlib.md5(url.encode(encoding='utf-8')).hexdigest()
            # if md5_name != 'f8c988df481e3046f1f90c81a1b98d90':
            #     continue
            if Path(f"{self.output}/{md5_name}.html").exists():
                logger.debug(f"{md5_name}.html exists, continue")
                continue
            new_url = "http://dx.doi.org/" + doi if doi else url
            logger.debug(f"{new_url} {md5_name}")
            driver.get(new_url)
            retry_num = 1
            try:
                while True:
                    time.sleep(10)
                    if url.split("=")[-1] == "Elsevier":
                        wait.until(ec.visibility_of_element_located((By.CSS_SELECTOR, '#gh-branding > svg')))
                        time.sleep((random.random() + 1) * 10)
                    elif "Wiley" in url.split("=")[-1]:
                        wait.until(ec.visibility_of_element_located((By.CSS_SELECTOR, '#mainLogo')))
                    elif "Thieme" in url.split("=")[-1]:
                        full_text = wait.until(
                            ec.element_to_be_clickable((By.CSS_SELECTOR, '#articleTabs > li:nth-child(2) > a')))
                        full_text.click()
                    elif "Pharmaceutical+Society+of+Japan" in url.split("=")[-1]:
                        full_text = wait.until(ec.element_to_be_clickable(
                            (By.CSS_SELECTOR, '#third-level-navigation > a.thirdlevel-active-btn')))
                        full_text.click()
                    elif url.split("=")[-1] in ["CHEMICAL+JOURNAL+OF+CHINESE+UNIVERSITIES-CHINESE",
                                                "Shaghai+Institute+of+Organic+Chemistry",
                                                "CHEMICAL+JOURNAL+OF+CHINESE+UNIVERSITIES-CHINESE"]:
                        full_text = wait.until(ec.element_to_be_clickable((By.CSS_SELECTOR,
                                                                           "#goTop > div.container.whitebg > "
                                                                           "div.abs-con > div > div > "
                                                                           "div.group.clearfix > div > "
                                                                           "div:nth-child(1) > span > a")))
                        full_text.click()
                    time.sleep((random.random() + 1) * 3)
                    try:
                        scroll_height = driver.execute_script("var q=document.body.scrollHeight; return(q)")
                    except JavascriptException:
                        logger.info(f"Loading {md5_name}.html failed, retry {retry_num}")
                        retry_num += 1
                        if retry_num >= 10:
                            break
                        else:
                            continue
                    else:
                        break
            except TimeoutException:
                logger.info(f"{md5_name}.html site error, continue")
                with open(f"{self.output}/{md5_name}.html", "w", encoding='utf-8') as f:
                    f.write(driver.page_source)
                continue

            time.sleep((random.random() + 1))

            # generate random scroll distance, total scroll_t times, sum is scrollHeight
            scroll_t = random.randint(20, 40)
            loop_height = np.random.dirichlet(np.ones(scroll_t)) * scroll_height
            for i, h in zip(range(scroll_t), loop_height):
                driver.execute_script(f"window.scrollBy(0,{h})")
                time.sleep(random.random()) if url.split("=")[-1] != "Elsevier" else time.sleep(random.random() + 1)
            driver.execute_script(f"window.scrollTo(0,{scroll_height})")
            time.sleep(random.random() * 3)

            # write to *.html files
            with open(f"{self.output}/{md5_name}.html", "w", encoding='utf-8') as f:
                f.write(driver.page_source)
            logger.info(f"{index:6d}: {md5_name}.html download successful")
            time.sleep(random.random() * 5)
            # exit()

        logger.info("All literatures download successful")
        driver.quit()

    def check_htmls(self):
        """
        Check the htmls crawled effectivity, maintain three files: exclude, checked and md5_name.json (from WOSParser)

        ----------------------------------------------------------------------
        md5_name.json  |   a mapping dict from html's name to its url
        exclude        |   the html in which is invalid but the inner error
        checked        |   the html is effective
        ----------------------------------------------------------------------
        """

        # load excluded
        with open(DataDir / "exclude", "r", encoding="utf-8") as f:
            exclude = [item.split(",")[0] for item in f.readlines()]

        # load checked
        with open(DataDir / "checked", "r", encoding="utf-8") as f:
            checked = [item.strip() for item in f.readlines()]

        # load md5_dict
        with open(DataDir / "md5_name.json", "r", encoding="utf-8") as f:
            md5_dict = json.load(f)

        # search which html failed
        literatures = Path(self.output).rglob("*.html")
        logger.info(f"Check download literature completeness")
        for literature in literatures:
            if literature.stem in checked:
                continue
            conclusion = []
            with open(literature, "r", encoding="utf-8") as f:
                html = f.read()
            try:
                doc = PyQuery(html)
                conclusion = re.findall(
                    "[Cc]onclusion|CONCLUSION|[Ss]ummary|[Oo]utlook|[Dd]iscussion|[Ee]xperiment|总结|结论|综上",
                    doc.text())
            except ParserError:
                pass

            if not len(conclusion) and literature.stem not in exclude:
                logger.info(f"{literature}: Search keyword failed, please check manually")
                logger.info(md5_dict[literature.stem.strip()])
            elif len(conclusion):
                with open(DataDir / "checked", "a+", encoding="utf-8") as f:
                    f.write(f"{literature.stem} \n")
                logger.info(f"{literature} store in database")


if __name__ == '__main__':
    crawler = PublisherCrawler(index_file=DataDir / "datafile.csv", output="../../literature")
    crawler.get_htmls()
    # crawler.check_htmls()
