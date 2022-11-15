from itertools import zip_longest

import logging
import numpy as np
from pyquery import PyQuery

from AICatalysis.common.file import HtmlIO

logger = logging.getLogger(__name__)

if __name__ == '__main__':
    files = ["0a27913b8f8012e7ef719ef175978842.html", "0abda187a9c883fcab7a0f45826ba61e.html",
             "0aed9192373c0cf784250adb59f7c9b6.html", "0b599ace8153b094978ce8de6d1b90a7.html",
             "0b900e98361640e82d938116ea3f9adc.html", "0b6079b37eb029d14a852b38d5ec4010.html",
             "0b3155901fc7c96509c7776641cc0964.html"]
    html_file = "../../literature/" + files[-1]
    html = HtmlIO.read(html_file)
    doc = PyQuery(html)
    tb = doc.find('div[class=article-table-content]')
    header = tb("header")
    if not len(tb):
        tb = doc.find('div.tableWrapper')
        header = tb("caption")

    tb_index = []
    for index, item in enumerate(header.items()):
        if "condition" in item.text().lower() or "phenylboronate" in item.text().lower():
            tb_index.append(index)

    header = header.filter(lambda i: i in tb_index).text()
    content = tb.filter(lambda i: i in tb_index)
    thead = content("thead")
    tbody = content("tbody")

    if len(thead("tr")) == 0:
        logger.warning("Can't find Table, Please check the html!!")
        exit(1)
    elif len(thead("tr")) == 1:
        thead_list = thead.text().splitlines()
    elif len(thead("tr")) == 2:  # solving many yields in two rows, e.g., yield (2a, 2b, 2c)
        row = []
        for tr_item in thead.items("tr"):
            if len(tr_item.find("th")):
                row.append([th_item.text() for index, th_item in enumerate(tr_item.items("th"))])
            elif len(tr_item.find("td")):
                row.append([th_item.text() for index, th_item in enumerate(tr_item.items("td"))])
        thead_list = []
        for item1, item2 in zip_longest(*row):
            if len(item2):
                if item1 is not None:
                    notation = item1
                    if len(item1):  # solving the figure row
                        thead_list.append(item1 + "-" + item2)
                    else:
                        thead_list.append(item2)
                else:
                    if len(notation):
                        thead_list.append(notation + "-" + item2)
                    else:
                        thead_list.append(item2)
            else:
                thead_list.append(item1)
    else:
        raise NotImplementedError

    tbody_list = np.array(tbody.text().splitlines()).reshape((-1, len(thead_list))).tolist()
    if len(tb.find("div.article-section__table-footnotes")):  # wily
        footnote = tb("div.article-section__table-footnotes").text()
    elif len(thead.parent().siblings("p")):  # www.thieme-connect.de
        footnote = thead.parent().siblings("p").text()

    # ----- print information ----- #
    print(header)
    print()

    print("+" + "-".center(25 * len(thead_list) + 2 * (len(thead_list) - 1), "-") + "+")
    for tr_item in thead_list:
        print("|" + f"{tr_item}".center(25) + "|", end="")
    print()
    print("|" + "-".center(25 * len(thead_list) + 2 * (len(thead_list) - 1), "-") + "|")
    for line in tbody_list:
        for tr_item in line:
            print("|" + f"{tr_item}".center(25) + "|", end="")
        print()
    print("+" + "-".center(25 * len(thead_list) + 2 * (len(thead_list) - 1), "-") + "+")

    print()
    print(footnote)
    pass
