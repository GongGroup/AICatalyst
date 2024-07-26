import abc
import copy
import logging
from itertools import zip_longest
from pathlib import Path

import numpy as np
import requests
from pyquery import PyQuery

from AICatalysis.common.constant import MD5Name, Table
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
        def tfchoose(item_):
            CODs = ['catalytic', 'catalysts', 'catalyzed', 'carbonylative', 'carbonylation', 'comparison', 'condition',
                    'control', 'development', 'effect', 'hydroformylation', 'important', 'influence', 'mechanistic',
                    'optimization', 'other', 'parameters', 'phenylboronate', 'poisoning', 'process', 'promoters',
                    'ratio', 'screen', 'selection', 'source', 'synthesis', 'various', 'with']
            for cod in CODs:
                if cod in item_.text().lower():
                    return True
            else:
                return False

        tb_index = []
        for index, item in enumerate(header.items()):
            if tfchoose(item):
                tb_index.append(index)
        return tb_index

    def _parse_table(self, tb=None, caption=None, th="th"):
        if not len(tb):
            logger.info("There's no table in literature, exit now!")
            exit(1)
        tb_index = BasePub.search_table(caption)  # obtain the index of the valid tables
        caption = caption.filter(lambda i: i in tb_index)  # type -> PyQuery
        table = tb.filter(lambda i: i in tb_index)  # type -> PyQuery

        if len(table):
            logger.info(f"Match {len(table)} tables, start parse~")
        else:
            logger.warning(f"No match table, please check!")

        def parse_single_table(stb_: PyQuery):
            """
            parse single table

            Args:
                stb_: single table prepare to parse (PyQuery)

            Returns:
                thead_list, tbody_list
            """
            thead = stb_("thead")  # type -> PyQuery
            tbody = stb_("tbody")  # type -> PyQuery

            # parse table-head
            if len(thead("tr")) == 0:
                thead_list_ = []
            elif len(thead("tr")) == 1:
                thead_list_ = thead.text().splitlines()
            elif len(thead("tr")) == 2:  # solving many yields in two rows, e.g., yield (2a, 2b, 2c)
                row = []
                for tr_item in thead.items("tr"):
                    row.append([th_item.text() for index, th_item in enumerate(tr_item.items(th))])

                if row[1][0]:
                    for _ in range(len(row[0]) - 1):
                        row[1].insert(0, '')

                two_row_list = []
                for item1, item2 in zip_longest(*row):
                    if item2 is not None and len(item2):
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
                        if item1 is not None:
                            two_row_list.append(item1)
                        else:
                            two_row_list.append("")
                thead_list_ = two_row_list
            elif len(thead("tr")) == 3:  # solving thead in three rows, e.g., 10.1002/ejoc.201201708
                thead_list_ = []
                tr_items = thead("tr").filter(
                    lambda i, this: not PyQuery(this).find("figure") and not PyQuery(this).find("img"))
                if len(tr_items) == 2:
                    tr_1 = tr_items.filter(lambda i: i == 0)
                    tr_2 = tr_items.filter(lambda i: i == 1)
                    for th_1, th_2 in zip(tr_1.items("th"), tr_2.items("th")):
                        th_1.remove("span")
                        th_2.remove("span")
                        thead_list_.append(th_1.text() + th_2.text())
                else:  # https://www.sciencedirect.com/science/article/pii/S0010854505000329?via%3Dihub
                    temp_thead_list_ = []
                    for trr in tr_items.items():
                        trr_list = []
                        for child in trr.children().items():
                            span = child.attr('colspan')
                            span = 1 if span is None else span
                            trr_list += [child.text()] * int(span)
                        temp_thead_list_.append(trr_list)
                    thead_list_ = [(i + "-" + j + "-" + k) for i, j, k in zip(*temp_thead_list_)]
            else:
                raise NotImplementedError

            # parse table-body
            if len(thead_list_):
                try:
                    tbody_list_ = np.array(tbody.text().splitlines()).reshape((-1, len(thead_list_))).tolist()
                except ValueError:  # '' in tbody, columns for each row is not match
                    tbody_list_ = []
                    row_span = 1
                    row_span_i, row_span_v = 0, None
                    for tr_pq in tbody.items("tr"):
                        tbody_row_list = []
                        for td_i, td_pq in enumerate(tr_pq.items("td")):
                            if row_span != 1 and td_i == row_span_i:
                                tbody_row_list.append(row_span_v)
                                row_span -= 1
                            if td_pq.attr("rowspan"):  # <td rowspan=3>...</td>
                                row_span = int(td_pq.attr("rowspan"))
                                row_span_i = td_i
                                row_span_v = td_pq.text()
                            tbody_row_list.append(td_pq.text())
                        tbody_list_.append(tbody_row_list)
            else:
                tbody_row_list = []
                for tr_pq in tbody.items("tr"):
                    tbody_row_list.append(tr_pq.text().split("\n"))
                tbody_list_ = tbody_row_list

            return thead_list_, tbody_list_  # type -> list, list

        thead_list, tbody_list = [], []
        for stb in table.items():
            stb_thead_list, stb_tbody_list = parse_single_table(stb)
            thead_list.append(stb_thead_list)
            tbody_list.append(stb_tbody_list)

        # parse table-caption and table-footnote
        footnote_text, caption_text = [], []
        for caption_pq in caption.items():  # type -> PyQuery
            caption_text.append(caption_pq.text())
            if len(caption_pq.parent()('.footnotes')):  # type -> PyQuery
                footnote_pq = caption_pq.parent()('.footnotes')
            else:
                footnote_pq = caption_pq.parent()('.legend').siblings('p')
            footnote_text.append(footnote_pq.text())

        assert len(caption_text) == len(thead_list) == len(tbody_list) == len(footnote_text)

        self._table = Table(caption_text, thead_list, tbody_list, footnote_text)  # type -> list, list, list, list

        return table  # type -> PyQuery

    @abc.abstractmethod
    def parse_table(self, printout=True, save=False, name="table.csv", url=None):
        if printout:
            self.print_table()
        if save:
            self.save_table(name=name, url=url)

    def print_table(self):
        for caption, thead_list, tbody_list, footnote in zip(*self._table):
            if not len(tbody_list):
                continue

            if not len(thead_list):
                thead_list = [''] * len(tbody_list[0])

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
                if not len(tbody):
                    continue
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
        self._table = Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)

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
            notation = ''
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

        table = self._parse_table(tb=tb, caption=header)

        # rewrite the parse table-footnote
        footnote_text = []
        for stb in table.items():
            tbody = stb("thead").parent()
            footnote_pq = tbody("tfoot")
            footnote_text.append(footnote_pq.text())
        self._table = Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)

        super(RSCPub, self).parse_table(**kargs)


class ACSPub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc.find('div[class=NLM_table-wrap]')
        header = tb("div[class=NLM_caption]")
        tb = self._parse_table(tb=tb, caption=header)

        # rewrite the parse table-footnote
        footnote_text = []
        for tb_pq in tb.items():
            if len(tb_pq("div[class=footnote]")):
                footnote_pq = tb_pq("div[class=footnote]")
            else:
                footnote_pq = tb_pq("div[class=NLM_table-wrap-foot]")
            footnote_text.append(footnote_pq.text())
        self._table = Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)

        super(ACSPub, self).parse_table(**kargs)


class ThiemePub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc.find('div.tableWrapper')
        header = tb("caption")

        table = self._parse_table(tb=tb, caption=header, th="td")
        # rewrite the parse table-footnote
        footnote_text = []
        for stb in table.items():
            thead = stb("thead")
            footnote_pq = thead.parent().siblings("p")
            footnote_text.append(footnote_pq.text())
        self._table = Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)

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
        tb = self.doc('div.table-wrap').parent('figure')
        caption = tb("figcaption")
        table = self._parse_table(tb=tb, caption=caption)

        # rewrite the parse table-footnote
        footnote_text = [stb("div.notes").text() for stb in table.items()]
        self._table = Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)
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
        atb = self.doc('div.c-article-table')
        url = atb('a')
        tb, caption = [], []
        for surl in url.items():
            href = surl.attr("href")
            real_href = 'https://link.springer.com' + href
            doc = PyQuery(requests.get(real_href).text)
            tb.append(doc('main'))
            caption.append(tb[-1]('h1'))
        table = self._parse_table(tb=PyQuery(tb), caption=PyQuery(caption))

        # rewrite the parse table-footnote
        footnote_text = [stb("div.c-article-table-footer").text() for stb in table.items()]
        self._table = Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)
        super(SpringerPub, self).parse_table(**kargs)


class KoreaPub(BasePub):

    def parse_table(self, **kargs):
        tb = self.doc.find('table').parent('div').parent('div')
        caption = tb(".captions")
        self._parse_table(tb=tb, caption=caption)
        super(KoreaPub, self).parse_table(**kargs)


class PharSocJapanPub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc.find('table').parent('div').parent('div')
        caption = PyQuery(copy.deepcopy(tb.copy()))  # remove the children node
        caption.find('div').remove()
        table = self._parse_table(tb=tb, caption=caption)

        # rewrite the parse table-footnote
        footnote_text = [stb("div:last-child").text() for stb in table.items()]
        self._table = Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)
        super(PharSocJapanPub, self).parse_table(**kargs)


class MDPIPub(BasePub):
    def parse_table(self, **kargs):
        if len(self.doc('div .html-caption')):
            tb = self.doc('div .html-table_show')
            caption = tb('div .html-caption')
            if not len(caption):
                caption = tb('caption')
            table = self._parse_table(tb=tb, caption=caption)

            # rewrite the parse table-footnote
            footnote_text = [stb('div .html-table_foot').text() for stb in table.items()]
            self._table = Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)
            super(MDPIPub, self).parse_table(**kargs)
        else:
            tb = self.doc('table')
            caption = tb('caption')
            self._parse_table(tb=tb, caption=caption)
            super(MDPIPub, self).parse_table(**kargs)


class BeiInsPub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc('table').parent('div').parent('figure')
        caption = tb('figcaption')
        table = self._parse_table(tb=tb, caption=caption)

        # rewrite the parse table-footnote
        footnote_text = [stb.next('p').text() for stb in table.items()]
        self._table = Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)
        super(BeiInsPub, self).parse_table(**kargs)


class NCBIPub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc('div.table-wrap')
        caption = tb('div.caption')
        table = self._parse_table(tb=tb, caption=caption)

        # rewrite the parse table-footnote
        footnote_text = [stb('div.tblwrap-foot').text() for stb in table.items()]
        self._table = Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)
        super(NCBIPub, self).parse_table(**kargs)


class IopSciPub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc('table').parent('div')
        caption = self.doc('table').prev('p')
        table = self._parse_table(tb=tb, caption=caption)

        # rewrite the parse table-footnote
        footnote_text = [stb('table').next('p').text() for stb in table.items()]
        self._table = Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)
        super(IopSciPub, self).parse_table(**kargs)


class HindawiPub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc('table').parent('div').parent('div')
        caption = tb("div.floats-partial-footer")
        self._parse_table(tb=tb, caption=caption)
        super(HindawiPub, self).parse_table(**kargs)


class NaturePub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc('div.c-article-table')
        caption = tb("figcaption")
        self._parse_table(tb=tb, caption=caption)
        super(NaturePub, self).parse_table(**kargs)


class CSJPub(BasePub):
    def parse_table(self, **kargs):
        tb = self.doc('div.NLM_table-wrap')
        caption = tb("div.caption")
        table = self._parse_table(tb=tb, caption=caption)

        # rewrite the parse table-footnote
        footnote_text = [stb('div.NLM_table-wrap-foot').text() for stb in table.items()]
        self._table = Table(self._table.caption, self._table.thead, self._table.tbody, footnote_text)
        super(CSJPub, self).parse_table(**kargs)


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
        'Pharmaceutical+Society+of+Japan': PharSocJapanPub,
        'MDPI+AG': MDPIPub,
        'Beilstein-Institut': BeiInsPub,
        'SAGE+Publications': SagePub,
        'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7452727': NCBIPub,
        'The+Electrochemical+Society': IopSciPub,
        'Hindawi+Limited': HindawiPub,
        'Springer+Science+and+Business+Media+LLC': NaturePub,
        'The+Chemical+Society+of+Japan': CSJPub,
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
    pass
