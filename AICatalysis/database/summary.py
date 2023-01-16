from pathlib import Path

from matplotlib import pyplot as plt

from AICatalysis.database.table_parser import HtmlTableParser


def get_pubs():
    literature_dir = "../../literature/"
    files = [file for file in Path(literature_dir).iterdir()]
    pubs = []
    for html_file in files:
        try:
            parser = HtmlTableParser(html_file)
        except KeyError:
            continue
        pubs.append(parser.pub)


def plot_pubs():
    # [('Elsevier', 435),
    #  ('Wiley+%28John+Wiley+%26+Sons%29', 293),
    #  ('American+Chemical+Society', 250),
    #  ('The+Royal+Society+of+Chemistry', 167),
    #  ('Springer-Verlag', 77),
    #  ('Georg+Thieme+Verlag+KG', 38),
    #  ('The+Chemical+Society+of+Japan', 33),
    #  ('MDPI+AG', 24),
    #  ('Informa+UK+%28Taylor+%26+Francis%29', 23)]
    plt.rc('font', family="Arial", weight="regular")
    data = [435, 293, 250, 167, 77, 38, 33, 24, 23, 260]
    plt.bar(['Elsevier', 'Wiley', 'ACS', 'RSC', 'Springer', 'Thieme', 'JCS', 'MDPI', 'Taylor', 'others'], data)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.show()


plot_pubs()
