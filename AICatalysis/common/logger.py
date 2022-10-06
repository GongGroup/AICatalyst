import logging
from pathlib import Path

from AICatalysis.common.constant import DATE, RESET, BOLD, YELLOW, WHITE, BLUE, RED, LogDir

COLORS = {
    'WARNING': YELLOW,
    'INFO': WHITE,
    'DEBUG': BLUE,
    'CRITICAL': YELLOW,
    'ERROR': RED
}


def formatter_message(message, use_color=True):
    if use_color:
        message = message.replace("$RESET", RESET).replace("$BOLD", BOLD)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")

    return message


class ColoredFormatter(logging.Formatter):

    def __init__(self, msg, use_color=False):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname

        if self.use_color and levelname in COLORS:
            levelname_color = COLORS[levelname] + levelname + RESET
            record.levelname = levelname_color

        return logging.Formatter.format(self, record)


class ColoredLogger(logging.Logger):
    FORMAT = "%(asctime)s ($BOLD%(filename)s$RESET:%(lineno)d) [$BOLD%(name)s$RESET][%(levelname)s] %(message)s"
    COLOR_FORMAT = formatter_message(FORMAT, True)
    FILE_FORMAT = formatter_message(FORMAT, False)

    def __init__(self, name):
        logging.Logger.__init__(self, name)
        file_formatter = ColoredFormatter(self.FILE_FORMAT, False)
        color_formatter = ColoredFormatter(self.COLOR_FORMAT, True)

        Path(LogDir).mkdir(exist_ok=True)
        fh = logging.FileHandler(f"{LogDir}/{DATE}.txt", encoding='utf-8')
        fh.setFormatter(file_formatter)

        ch = logging.StreamHandler()
        ch.setFormatter(color_formatter)

        self.addHandler(fh)
        self.addHandler(ch)


def init_root_logger(name=None, level=logging.INFO):
    logging.setLoggerClass(ColoredLogger)
    root_logger = logging.getLogger(name)
    root_logger.setLevel(level=level)
