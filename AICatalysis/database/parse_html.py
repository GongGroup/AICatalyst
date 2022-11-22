import abc
import logging
from itertools import zip_longest
from pathlib import Path

import numpy as np
from pyquery import PyQuery

from AICatalysis.common.constant import MD5Name, _Table
from AICatalysis.common.file import HtmlIO, JsonIO
from AICatalysis.common.logger import init_root_logger

init_root_logger("AICatalysis")
logger = logging.getLogger(__name__)


class BasePub(metaclass=abc.ABCMeta):
    def __init__(self, file):
        self._html = HtmlIO.read(file)
        self.doc = PyQuery(self._html)

        self._table = None

    @abc.abstractmethod
    def parse_table(self):
        pass

    def _parse_table(self, tb=None, caption=None, th="th"):
        tb_index = BasePub.search_table(caption)  # obtain the index of the valid tables
        caption = caption.filter(lambda i: i in tb_index)  # type -> PyQuery
        content = tb.filter(lambda i: i in tb_index)  # type -> PyQuery
        thead = content("thead")  # type -> PyQuery
        tbody = content("tbody")  # type -> PyQuery
        logger.info(f"Match {len(thead)} tables, start parse~")

        # parse table-head
        thead_list = []
        for thead_pq in thead.items():  # type -> PyQuery
            if len(thead_pq("tr")) == 0:
                logger.warning("Can't find Table, Please check the html!!")
                exit(1)
            elif len(thead_pq("tr")) == 1:
                thead_list.append(thead_pq.text().splitlines())
            elif len(thead_pq("tr")) == 2:  # solving many yields in two rows, e.g., yield (2a, 2b, 2c)
                row = []
                for tr_item in thead_pq.items("tr"):
                    row.append([th_item.text() for index, th_item in enumerate(tr_item.items(th))])

                if row[1][0]:
                    for _ in range(len(row[0]) - 1):
                        row[1].insert(0, '')

                two_row_list = []
                for item1, item2 in zip_longest(*row):
                    if len(item2):
                        if item1 is not None:
                            notation = item1
                            if len(item1):  # solving the figure row
                                two_row_list.append(item1 + "-" + item2)
                            else:
                                two_row_list.append(item2)
                        else:
                            if len(notation):
                                two_row_list.append(notation + "-" + item2)
                            else:
                                two_row_list.append(item2)
                    else:
                        two_row_list.append(item1)
                thead_list.append(two_row_list)
            else:
                raise NotImplementedError

        # parse table-content
        tbody_list = []
        for tbody_pq, thead_pq in zip(tbody.items(), thead_list):  # type -> PyQuery, list
            tbody_pq_list = np.array(tbody_pq.text().splitlines()).reshape((-1, len(thead_pq))).tolist()
            tbody_list.append(tbody_pq_list)

        # parse table-caption and table-footnote
        footnote_text, caption_text = [], []
        for caption_pq in caption.items():  # type -> PyQuery
            caption_text.append(caption_pq.text())
            footnote_pq = caption_pq.parent()('.footnotes')  # type -> PyQuery
            footnote_text.append(footnote_pq.text())

        self._table = _Table(caption_text, thead_list, tbody_list, footnote_text)  # type -> list, list, list, list

        return thead  # type -> PyQuery

    @staticmethod
    def search_table(header):
        def tfchoose(item):
            CODs = ['condition', 'phenylboronate', 'screen', 'effect']
            for cod in CODs:
                if cod in item.text().lower():
                    return True
            else:
                return False

        tb_index = []
        for index, item in enumerate(header.items()):
            if tfchoose(item):
                tb_index.append(index)
        return tb_index

    def print_table(self):

        for caption, thead_list, tbody_list, footnote in zip(*self._table):
            print(caption)
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
            print()


class WileyPub(BasePub):
    def parse_table(self):
        tb = self.doc.find('div[class=article-table-content]')
        header = tb("header")

        tb_index = WileyPub.search_table(header)
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
                row.append([th_item.text() for index, th_item in enumerate(tr_item.items("th"))])

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
        footnote = tb("div.article-section__table-footnotes").text()

        self._parse_table = header, thead_list, tbody_list, footnote
        self.print_table()


class SJOCPub(BasePub):
    def parse_table(self):
        tb = self.doc.find('div[class=article-table-content]')
        header = tb("header")

        tb_index = WileyPub.search_table(header)
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
                row.append([th_item.text() for index, th_item in enumerate(tr_item.items("th"))])

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
        footnote = tb("div.article-section__table-footnotes").text()

        self._parse_table = header, thead_list, tbody_list, footnote
        self.print_table()


class RSCPub(BasePub):
    def parse_table(self):
        tb = self.doc.find('div[class=article-table-content]')
        header = tb("header")

        tb_index = WileyPub.search_table(header)
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
                row.append([th_item.text() for index, th_item in enumerate(tr_item.items("th"))])

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
        footnote = tb("div.article-section__table-footnotes").text()

        self._parse_table = header, thead_list, tbody_list, footnote
        self.print_table()


class ACSPub(BasePub):
    def parse_table(self):
        tb = self.doc.find('div[class=article-table-content]')
        header = tb("header")

        tb_index = WileyPub.search_table(header)
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
                row.append([th_item.text() for index, th_item in enumerate(tr_item.items("th"))])

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
        footnote = tb("div.article-section__table-footnotes").text()

        self._parse_table = header, thead_list, tbody_list, footnote
        self.print_table()


class ThiemePub(BasePub):
    def parse_table(self):
        tb = self.doc.find('div.tableWrapper')
        header = tb("caption")

        thead = self._parse_table(tb=tb, caption=header, th="td")
        # rewrite the parse table-footnote
        footnote_text = []
        for thead_pq in thead.items():
            footnote_pq = thead_pq.parent().siblings("p")
            footnote_text.append(footnote_pq.text())
        self._table = _Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)

        self.print_table()


class TaylorPub(BasePub):
    def parse_table(self):
        tb = self.doc.find('div.tableWrapper')
        header = tb("caption")

        tb_index = BasePub.search_table(header)
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
        footnote = thead.parent().siblings("p").text()

        self._parse_table = header, thead_list, tbody_list, footnote
        self.print_table()


class ElsevierPub(BasePub):
    def parse_table(self):
        tb = self.doc.find('table').parent('div').parent('div')
        caption = tb(".captions")
        self._parse_table(tb=tb, caption=caption)
        self.print_table()


class HtmlTableParser(object):
    Allocator = {
        "American+Chemical+Society": ACSPub,
        "Elsevier": ElsevierPub,
        "The+Royal+Society+of+Chemistry": RSCPub,
        "Shaghai+Institute+of+Organic+Chemistry": SJOCPub,
        "Informa+UK+%28Taylor+%26+Francis%29": TaylorPub,
        "Georg+Thieme+Verlag+KG": ThiemePub,
        "Wiley+%28John+Wiley+%26+Sons%29": WileyPub,
    }

    def __init__(self, file):
        self.file = file
        self.name = Path(self.file).stem
        self.pub = JsonIO.read(MD5Name)[self.name].split("=")[-1]

    def parse(self):
        return HtmlTableParser.Allocator[self.pub](self.file).parse_table()


if __name__ == '__main__':
    files = ["0a27913b8f8012e7ef719ef175978842.html", "0abda187a9c883fcab7a0f45826ba61e.html",
             "0aed9192373c0cf784250adb59f7c9b6.html", "0b599ace8153b094978ce8de6d1b90a7.html",
             "0b900e98361640e82d938116ea3f9adc.html", "0b6079b37eb029d14a852b38d5ec4010.html",
             "0b3155901fc7c96509c7776641cc0964.html", "0c919f91371e15e113a3cc772518c02b.html",
             "0c9734f494e6d610f9e8d44879db1e70.html"]
    html_file = "../../literature/" + files[6]

    parser = HtmlTableParser(html_file)
    parser.parse()

    pass
