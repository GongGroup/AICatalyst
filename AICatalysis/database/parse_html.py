import abc
import logging
from itertools import zip_longest
from pathlib import Path

import numpy as np
from pyquery import PyQuery

from AICatalysis.common.constant import MD5Name, _Table
from AICatalysis.common.file import HtmlIO, JsonIO
from AICatalysis.common.logger import init_root_logger

init_root_logger()
logger = logging.getLogger(__name__)


class BasePub(metaclass=abc.ABCMeta):
    def __init__(self, file):
        self._html = HtmlIO.read(file)
        self.doc = PyQuery(self._html)

        self._table = None

    @staticmethod
    def search_table(header):
        def tfchoose(item):
            CODs = ['condition', 'phenylboronate', 'screen', 'effect', 'optimization', 'control', 'various',
                    'with and without']
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

    def _parse_table(self, tb=None, caption=None, th="th"):
        tb_index = BasePub.search_table(caption)  # obtain the index of the valid tables
        caption = caption.filter(lambda i: i in tb_index)  # type -> PyQuery
        content = tb.filter(lambda i: i in tb_index)  # type -> PyQuery
        thead = content("thead")  # type -> PyQuery
        tbody = content("tbody")  # type -> PyQuery
        logger.info(f"Match {len(thead)} tables, start parse~") if len(thead) \
            else logger.warning(f"Search table failed")

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
            try:
                tbody_pq_list = np.array(tbody_pq.text().splitlines()).reshape((-1, len(thead_pq))).tolist()
            except ValueError:  # '' in tbody, columns for each row is not match
                tbody_row_list = []
                for td_pq in tbody_pq.items("td"):
                    tbody_row_list.append(td_pq.text())
                tbody_pq_list = np.array(tbody_row_list).reshape((-1, len(thead_pq))).tolist()
            tbody_list.append(tbody_pq_list)

        # parse table-caption and table-footnote
        footnote_text, caption_text = [], []
        for caption_pq in caption.items():  # type -> PyQuery
            caption_text.append(caption_pq.text())
            footnote_pq = caption_pq.parent()('.footnotes')  # type -> PyQuery
            footnote_text.append(footnote_pq.text())

        self._table = _Table(caption_text, thead_list, tbody_list, footnote_text)  # type -> list, list, list, list

        return thead  # type -> PyQuery

    @abc.abstractmethod
    def parse_table(self, print=True, save=False, name="table.csv", url=None):
        if print:
            self.print_table()
        if save:
            self.save_table(name=name, url=url)

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

    def save_table(self, name="table.csv", url=None):

        if not len(self._table[2]):
            logger.warning(f"There is no information in parsed table!!")
            return

        with open(name, "w", encoding="utf-8") as f:
            logger.info(f"Store table in `{name}`")
            for caption, thead, tbody, footnote in zip(*self._table):
                f.write(caption + "\n")
                f.write(",".join(thead) + "\n")
                for item in tbody:
                    f.write(",".join(item) + "\n")
                f.write(footnote + "\n")
                f.write(url + "\n")
                f.write("\n")


class WileyPub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc.find('div[class=article-table-content]')
        header = tb("header")

        self._parse_table(tb=tb, caption=header)

        # rewrite the parse table-footnote
        footnote_text = []
        for tb_pq in tb.items():
            footnote_pq = tb_pq("div.article-section__table-footnotes")
            footnote_text.append(footnote_pq.text())
        self._table = _Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)

        super(WileyPub, self).parse_table(**kargs)


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
    def parse_table(self, **kargs):
        tb = self.doc.find('figure[class=pnl--table]')
        header = tb("figcaption")

        thead = self._parse_table(tb=tb, caption=header)

        # rewrite the parse table-footnote
        footnote_text = []
        for tb_pq in thead.parent().items():
            footnote_pq = tb_pq("tfoot")
            footnote_text.append(footnote_pq.text())
        self._table = _Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)

        super(RSCPub, self).parse_table(**kargs)


class ACSPub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc.find('div[class=NLM_table-wrap]')
        header = tb("div[class=NLM_caption]")

        self._parse_table(tb=tb, caption=header)

        # rewrite the parse table-footnote
        footnote_text = []
        for tb_pq in tb.items():
            footnote_pq = tb_pq("div[class=footnote]")
            footnote_text.append(footnote_pq.text())
        self._table = _Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)

        super(ACSPub, self).parse_table(**kargs)


class ThiemePub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc.find('div.tableWrapper')
        header = tb("caption")

        thead = self._parse_table(tb=tb, caption=header, th="td")
        # rewrite the parse table-footnote
        footnote_text = []
        for thead_pq in thead.items():
            footnote_pq = thead_pq.parent().siblings("p")
            footnote_text.append(footnote_pq.text())
        self._table = _Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)

        super(ThiemePub, self).parse_table(**kargs)


class TaylorPub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc.find('div.tableWrapper')
        header = tb("caption")

        self._parse_table(tb=tb, caption=header)
        super(TaylorPub, self).parse_table(**kargs)


class ElsevierPub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc.find('table').parent('div').parent('div')
        caption = tb(".captions")
        self._parse_table(tb=tb, caption=caption)
        super(ElsevierPub, self).parse_table(**kargs)


class PlosPub(BasePub):

    def parse_table(self, **kargs):
        tb = self.doc.find('table').parent('div').parent('div')
        caption = tb(".captions")
        self._parse_table(tb=tb, caption=caption)
        super(PlosPub, self).parse_table(**kargs)


class SagePub(BasePub):

    def parse_table(self, **kargs):
        tb = self.doc.find('table').parent('div').parent('div')
        caption = tb(".captions")
        self._parse_table(tb=tb, caption=caption)
        super(SagePub, self).parse_table(**kargs)


class FunpecrpPub(BasePub):

    def parse_table(self, **kargs):
        tb = self.doc.find('table').parent('div').parent('div')
        caption = tb(".captions")
        self._parse_table(tb=tb, caption=caption)
        super(FunpecrpPub, self).parse_table(**kargs)


class EurekaselectPub(BasePub):

    def parse_table(self, **kargs):
        tb = self.doc.find('table').parent('div').parent('div')
        caption = tb(".captions")
        self._parse_table(tb=tb, caption=caption)
        super(EurekaselectPub, self).parse_table(**kargs)


class SpringerPub(BasePub):

    def parse_table(self, **kargs):
        tb = self.doc.find('table').parent('div').parent('div')
        caption = tb(".captions")
        self._parse_table(tb=tb, caption=caption)
        super(SpringerPub, self).parse_table(**kargs)


class KoreaPub(BasePub):

    def parse_table(self, **kargs):
        tb = self.doc.find('table').parent('div').parent('div')
        caption = tb(".captions")
        self._parse_table(tb=tb, caption=caption)
        super(KoreaPub, self).parse_table(**kargs)


class HtmlTableParser(object):
    Allocator = {
        "American+Chemical+Society": ACSPub,
        "Elsevier": ElsevierPub,
        "The+Royal+Society+of+Chemistry": RSCPub,
        "Shaghai+Institute+of+Organic+Chemistry": SJOCPub,
        "Informa+UK+%28Taylor+%26+Francis%29": TaylorPub,
        "Georg+Thieme+Verlag+KG": ThiemePub,
        "Wiley+%28John+Wiley+%26+Sons%29": WileyPub,
        "Public+Library+of+Science": PlosPub,
        "Science+Reviews+2000+LTD.": SagePub,
        "Genetics+and+Molecular+Research": FunpecrpPub,
        'Bentham+Science': EurekaselectPub,
        'Springer-Verlag': SpringerPub,
        'Pleiades+Publishing': SpringerPub,
        'Korean+Chemical+Society': KoreaPub,
    }

    def __init__(self, file):
        self.file = file
        self.name = Path(self.file).stem
        self.url = JsonIO.read(MD5Name)[self.name]
        self.pub = self.url.split("=")[-1]

    def parse(self, **kargs):
        logger.info(f"Start parse `{self.url}`")
        return HtmlTableParser.Allocator[self.pub](self.file).parse_table(**kargs)


if __name__ == '__main__':
    literature_dir = "../../literature/"
    files = [file for file in Path(literature_dir).iterdir()]
    html_file = files[41]

    parser = HtmlTableParser(html_file)
    parser.parse(save=True, name=f"{parser.name}.csv", url=parser.url)

    pass
