import numpy as np
from pyquery import PyQuery

from AICatalysis.common.file import HtmlIO

if __name__ == '__main__':
    html_file = "../../literature/0a27913b8f8012e7ef719ef175978842.html"
    html = HtmlIO.read(html_file)
    doc = PyQuery(html)
    tb1 = doc.find("#tbl1")
    header = tb1("header").text()

    content = tb1("div.article-table-content-wrapper")
    thead = content("thead")
    tbody = content("tbody")
    thead_list = thead.text().splitlines()
    tbody_list = np.array(tbody.text().splitlines()).reshape((-1, len(thead_list))).tolist()

    footnote = tb1("div.article-section__table-footnotes").text()

    print(header)
    print()
    print("+" + "-".center(25 * len(thead_list) + 2 * (len(thead_list) - 1), "-") + "+")
    for item in thead_list:
        print("|" + f"{item}".center(25) + "|", end="")
    print()
    print("|" + "-".center(25 * len(thead_list) + 2 * (len(thead_list) - 1), "-") + "|")
    for line in tbody_list:
        for item in line:
            print("|" + f"{item}".center(25) + "|", end="")
        print()
    print("+" + "-".center(25 * len(thead_list) + 2 * (len(thead_list) - 1), "-") + "+")

    print()
    print(footnote)
    pass
