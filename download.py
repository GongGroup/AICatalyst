import hashlib
import random
import time
from pathlib import Path

import numpy as np
from selenium import webdriver
import undetected_chromedriver as uc
from selenium.common import JavascriptException
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait

from logger import logger

if __name__ == '__main__':

    with open("datafile.1.csv", encoding='utf-8') as f:
        content = f.readlines()

    options = uc.ChromeOptions()
    options.add_argument(r"user-data-dir=C:\Users\hui_zhou\AppData\Local\Google\Chrome\User Data")
    options.page_load_strategy = 'none'
    driver = uc.Chrome(options=options)
    driver.set_window_position(1200, 10)
    wait = WebDriverWait(driver, 60)

    for line in content[:45]:
        url, doi = line.strip().split(",")[-2:]
        md5_name = hashlib.md5(url.encode(encoding='utf-8')).hexdigest()
        if Path(f"literature/{md5_name}.html").exists():
            logger.info(f"{md5_name}.html exists, continue")
            continue
        new_url = "http://dx.doi.org/" + doi if doi else url

        driver.get(new_url)
        while True:
            count = 1
            time.sleep(10)
            if url.split("=")[-1] == "Elsevier":
                time.sleep(random.random() * 30)
                wait.until(EC.visibility_of_element_located((By.CSS_SELECTOR, '#gh-branding > svg')))
                time.sleep(random.random() * 3)
            try:
                scrollHeight = driver.execute_script("var q=document.body.scrollHeight; return(q)")
            except JavascriptException:
                logger.info(f"Loading html failed, retry {count}")
                count += 1
                continue
            else:
                break

        # generate random scroll distance, total scroll_t times, sum is scrollHeight
        scroll_t = random.randint(20, 40)
        loop_height = np.random.dirichlet(np.ones(scroll_t)) * scrollHeight
        for i, h in zip(range(scroll_t), loop_height):
            driver.execute_script(f"window.scrollBy(0,{h})")
            time.sleep(random.random())
        driver.execute_script(f"window.scrollTo(0,{scrollHeight})")
        time.sleep(random.random() * 3)

        # write to *.html files
        with open(f"literature/{md5_name}.html", "w", encoding='utf-8') as f:
            f.write(driver.page_source)
        logger.info(f"{md5_name}.html download successful")
        time.sleep(random.random() * 5)

    logger.info("All literatures download successful")
    driver.quit()
