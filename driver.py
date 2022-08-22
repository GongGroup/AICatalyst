from selenium import webdriver
import undetected_chromedriver as uc
from selenium.webdriver.support.wait import WebDriverWait


class ChromeDriver(object):
    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super().__new__(cls)

        return cls._instance

    def __init__(self, timeout=120, options=True):
        self.timeout = timeout
        self.options = options
        self._driver = None

    @property
    def driver(self):
        if self._driver is None:
            if self.options:
                options = uc.ChromeOptions()
                options.add_argument(r"user-data-dir=C:\Users\hui_zhou\AppData\Local\Google\Chrome\User Data")
                options.page_load_strategy = 'none'
                self._driver = uc.Chrome(options=options)
            else:
                self._driver = uc.Chrome()
            self._driver.set_window_position(1200, 10)

        return self._driver

    @property
    def wait(self):
        return WebDriverWait(self.driver, self.timeout)
