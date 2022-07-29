import json
import random
import time

import numpy as np
from pyquery import PyQuery as pq
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait

from logger import logger

with open("config.json", "r") as f:
    config = json.load(f)

# load config setting
root_url = config['WOS_root']
start_page = config['start_page']
end_page = config['end_page']
keyword = config['keyword']

# start driver
driver = webdriver.Chrome()
driver.set_window_position(1200, 10)
wait = WebDriverWait(driver, 60)
driver.get(root_url)

logger.info("prepare to accept cookie")
cookie_accept = wait.until(EC.visibility_of_element_located((By.CSS_SELECTOR, '#onetrust-accept-btn-handler')))
driver.execute_script("arguments[0].click();", cookie_accept)

logger.info("prepare to input search")
search_input = wait.until(EC.visibility_of_element_located((By.XPATH, '//*[@id="mat-input-0"]')))
search_input.click()
search_input.send_keys(keyword)

logger.info("prepare to submit")
submit_button = wait.until(EC.element_to_be_clickable((By.XPATH, '//*[@id="snSearchType"]/div[3]/button[2]/span[1]')))
submit_button.click()
time.sleep(random.random() * 9)

# obtain total pages
cur_url = driver.current_url
time.sleep(random.random())
doc = pq(driver.page_source)
tot_pages = int(doc("span[class=end-page\ ng-star-inserted]").text().split()[0])
logger.info(f"Total Pages: {tot_pages}")

if start_page > tot_pages:
    raise ValueError(f"start_page = {start_page} over total_pages({tot_pages})")

if end_page > tot_pages:
    logger.warning(f"end_page = {end_page} over total_pages({tot_pages}), reset end_page = {tot_pages}")
    end_page = tot_pages

cur_page = start_page

# start downloading
while cur_page <= end_page:

    # set random sleep
    sleep_t1 = random.random() + 1
    sleep_t3 = sleep_t1 * 3
    sleep_t5 = sleep_t1 * 5
    sleep_t9 = sleep_t1 * 9

    logger.info(f"Page {cur_page} downloading...")
    cur_url = "/".join(cur_url.split("/")[:-1]) + f"/{cur_page}"
    cur_html = driver.get(cur_url)
    time.sleep(sleep_t3)

    # obtain the scrollHeight by js
    scrollHeight = driver.execute_script("var q=document.body.scrollHeight; return(q)")

    # generate random scroll distance, total scroll_t times, sum is scrollHeight
    scroll_t = random.randint(20, 40)
    loop_height = np.random.dirichlet(np.ones(scroll_t)) * scrollHeight
    for i, h in zip(range(scroll_t), loop_height):
        driver.execute_script(f"window.scrollBy(0,{h})")
        time.sleep(sleep_t1)
    driver.execute_script(f"window.scrollTo(0,{scrollHeight})")

    # write to *.html files
    with open(f"htmls/page_{cur_page}.html", "w", encoding='utf-8') as f:
        f.write(driver.page_source)
    cur_page += 1

    if cur_page % 5 == 0:
        time.sleep(sleep_t5)

    if cur_page % 20 == 0:
        time.sleep(sleep_t9)

logger.info("All Pages download successful")

driver.close()
