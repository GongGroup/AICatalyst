from itertools import zip_longest

import numpy as np
from pyquery import PyQuery

from AICatalysis.common.file import HtmlIO

if __name__ == '__main__':
    html_file = "../../literature/0abda187a9c883fcab7a0f45826ba61e.html"
    html = HtmlIO.read(html_file)
    doc = PyQuery(html)
    tb1 = doc.find("#tbl1")
    header = tb1("header").text()

    content = tb1("div.article-table-content-wrapper")
    thead = content("thead")
    tbody = content("tbody")

    if len(thead("tr")) == 1:
        thead_list = thead.text().splitlines()
    elif len(thead("tr")) == 2:  # solving many yields in two rows, e.g., yield (2a, 2b, 2c)
        row = []
        for tr_item in thead.items("tr"):
            row.append([th_item.text() for index, th_item in enumerate(tr_item.items("th"))])
        thead_list = []
        for item1, item2 in zip_longest(*row):
            if len(item2):
                if item1 is not None:
                    thead_list.append(item1 + "-" + item2)
                    notation = item1
                else:
                    thead_list.append(notation + "-" + item2)
            else:
                thead_list.append(item1)
    else:
        raise NotImplementedError

    tbody_list = np.array(tbody.text().splitlines()).reshape((-1, len(thead_list))).tolist()
    footnote = tb1("div.article-section__table-footnotes").text()

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
