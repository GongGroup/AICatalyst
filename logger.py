import os
import time
import inspect
import logging

CurrentDir = os.path.dirname(os.path.abspath(os.path.realpath(inspect.getfile(inspect.currentframe()))))
RootDir = CurrentDir
LogDir = os.path.join(CurrentDir, f"{RootDir}/logs")

if not os.path.exists(LogDir):
    os.mkdir(LogDir)

Date = time.strftime("%Y-%m-%d", time.localtime())
BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

COLORS = {
    'WARNING': YELLOW,
    'INFO': WHITE,
    'DEBUG': BLUE,
    'CRITICAL': YELLOW,
    'ERROR': RED
}


def formatter_message(message, use_color=True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
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
            levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
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

        fh = logging.FileHandler(f"{LogDir}/{Date}.txt")
        fh.setFormatter(file_formatter)

        ch = logging.StreamHandler()
        ch.setFormatter(color_formatter)

        self.addHandler(fh)
        self.addHandler(ch)


logging.setLoggerClass(ColoredLogger)
logger = logging.getLogger("__main__")
logger.setLevel(logging.INFO)
