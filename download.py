import hashlib
import random
import time
from pathlib import Path

import numpy as np
from selenium import webdriver

from logger import logger

with open("datafile.1.csv", encoding='utf-8') as f:
    content = f.readlines()

options = webdriver.ChromeOptions()
options.add_argument(r"user-data-dir=C:\Users\hui_zhou\AppData\Local\Google\Chrome\User Data")
options.add_argument("blink-settings=imagesEnabled=false")

driver = webdriver.Chrome(options=options)
driver.set_window_position(1200, 10)

for line in content[:5]:
    url, doi = line.strip().split(",")[-2:]
    md5_name = hashlib.md5(url.encode(encoding='utf-8')).hexdigest()
    if Path(f"literature/{md5_name}.html").exists():
        logger.info(f"{md5_name}.html exists, continue")
        continue

    new_url = "http://dx.doi.org/" + doi if doi else doi
    driver.get(new_url)
    time.sleep(random.random() * 3)

    scrollHeight = driver.execute_script("var q=document.body.scrollHeight; return(q)")

    # generate random scroll distance, total scroll_t times, sum is scrollHeight
    scroll_t = random.randint(20, 40)
    loop_height = np.random.dirichlet(np.ones(scroll_t)) * scrollHeight
    for i, h in zip(range(scroll_t), loop_height):
        driver.execute_script(f"window.scrollBy(0,{h})")
        time.sleep(random.random())
    driver.execute_script(f"window.scrollTo(0,{scrollHeight})")

    # write to *.html files
    with open(f"literature/{md5_name}.html", "w", encoding='utf-8') as f:
        f.write(driver.page_source)
    logger.info(f"{md5_name}.html download successful")
    time.sleep(random.random() * 5)

logger.info("All literatures download successful")
driver.quit()
