import hashlib
import logging
import random
import time
from pathlib import Path

import numpy as np
from selenium import webdriver
import undetected_chromedriver as uc
from selenium.common import JavascriptException
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait

from logger import logger

logger.setLevel(logging.INFO)
if __name__ == '__main__':

    with open("datafile.1.csv", encoding='utf-8') as f:
        content = f.readlines()

    options = uc.ChromeOptions()
    options.add_argument(r"user-data-dir=C:\Users\hui_zhou\AppData\Local\Google\Chrome\User Data")
    options.page_load_strategy = 'none'
    driver = uc.Chrome(options=options)
    driver.set_window_position(1200, 10)
    wait = WebDriverWait(driver, 120)
    scrollHeight = 1000

    for line in content[:800]:
        url, doi = line.strip().split(",")[-2:]
        md5_name = hashlib.md5(url.encode(encoding='utf-8')).hexdigest()
        if Path(f"literature/{md5_name}.html").exists():
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
                    wait.until(EC.visibility_of_element_located((By.CSS_SELECTOR, '#gh-branding > svg')))
                    time.sleep((random.random() + 1) * 10)
                elif "Wiley" in url.split("=")[-1]:
                    wait.until(EC.visibility_of_element_located((By.CSS_SELECTOR, '#mainLogo')))
                elif "Thieme" in url.split("=")[-1]:
                    full_text = wait.until(
                        EC.element_to_be_clickable((By.CSS_SELECTOR, '#articleTabs > li:nth-child(2) > a')))
                    full_text.click()
                elif url.split("=")[-1] in ["CHEMICAL+JOURNAL+OF+CHINESE+UNIVERSITIES-CHINESE",
                                            "Shaghai+Institute+of+Organic+Chemistry",
                                            "CHEMICAL+JOURNAL+OF+CHINESE+UNIVERSITIES-CHINESE"]:
                    full_text = wait.until(EC.element_to_be_clickable((By.CSS_SELECTOR,
                                                                       "#goTop > div.container.whitebg > div.abs-con > div "
                                                                       "> div > div.group.clearfix > div > div:nth-child(1) "
                                                                       "> span > a")))
                    full_text.click()
                time.sleep((random.random() + 1) * 3)
                try:
                    scrollHeight = driver.execute_script("var q=document.body.scrollHeight; return(q)")
                except JavascriptException:
                    logger.info(f"Loading {md5_name}.html failed, retry {retry_num}")
                    retry_num += 1
                    continue
                else:
                    break
        except TimeoutException:
            logger.info(f"{md5_name}.html site error, continue")
            with open(f"literature/{md5_name}.html", "w", encoding='utf-8') as f:
                f.write(driver.page_source)
            continue

        time.sleep((random.random() + 1))

        # generate random scroll distance, total scroll_t times, sum is scrollHeight
        scroll_t = random.randint(20, 40)
        loop_height = np.random.dirichlet(np.ones(scroll_t)) * scrollHeight
        for i, h in zip(range(scroll_t), loop_height):
            driver.execute_script(f"window.scrollBy(0,{h})")
            time.sleep(random.random()) if url.split("=")[-1] != "Elsevier" else time.sleep(random.random() + 1)
        driver.execute_script(f"window.scrollTo(0,{scrollHeight})")
        time.sleep(random.random() * 3)

        # write to *.html files
        with open(f"literature/{md5_name}.html", "w", encoding='utf-8') as f:
            f.write(driver.page_source)
        logger.info(f"{md5_name}.html download successful")
        time.sleep(random.random() * 5)

    logger.info("All literatures download successful")
    driver.quit()
