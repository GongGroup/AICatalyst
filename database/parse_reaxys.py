#input file
from common.constant import FReaxysXML, FRXconfig, FRXDconfig
#output file
from common.constant import FReaxys

from common.fio import JsonIO
from common.logger import logger
from bs4 import BeautifulSoup

class ParseReaxys:
    def __init__(self):
        self.rx_config = JsonIO.read(FRXconfig)
        self.rxd_config = JsonIO.read(FRXDconfig)
        self.data_type = []
        self.data_reaction = []
        self._initialize()

    def _initialize(self):
        with open(FReaxysXML, 'r') as f:
            self.reaxys_xml = f.read()

    def forward(self):
        logger.info('start to parse reaxys_xml')
        reaction_etree = BeautifulSoup(self.reaxys_xml, 'xml')
        reaction_type_list = reaction_etree.find_all(name='reaction')
        for index, reaction_type in enumerate(reaction_type_list):

            self._parse_rx(reaction_type)

            rxd_list = reaction_type.find_all('RXD')
            for cur_RXD in rxd_list:
                self._parse_rxd(cur_RXD)
            logger.info(f"reaction {index + 1} has been parsed")
        JsonIO.write(self.data_reaction, FReaxys)
        logger.info("Successfully parse reaxys_xml")

    def list2dict(func):
        def decorate(self, *args):
            if func.__name__ == '_parse_rx':
                key_list, value_list = func(self, *args)
                dict_type = dict(zip(key_list, value_list))
                self.data_type.append(dict_type)
            elif func.__name__ == '_parse_rxd':
                key_list, value_list = func(self, *args)
                dict_react = dict(zip(key_list, value_list))
                dict_react.update(self.data_type[-1])
                self.data_reaction.append(dict_react)
        return decorate

    @list2dict
    def _parse_rx(self, reaction_type):
        key_list = []
        value_list = []

        for event, elem in self.rx_config.items():
            key_list.append(event)

            elem_text = reaction_type.find_all(elem)
            elem_text_list = []
            for i in range(len(elem_text)):
                elem_text_list.append(elem_text[i].text)
            value_list.append(elem_text_list)
        return key_list, value_list

    @list2dict
    def _parse_rxd(self, cur_RXD):
        key_list = []
        value_list = []

        for event, elem in self.rxd_config.items():
            key_list.append(event)

            elem_value_list = []
            elem_list = cur_RXD.find_all(elem)
            for i in range(len(elem_list)):
                elem_value_list.append(elem_list[i].text)
            value_list.append(elem_value_list)
        return key_list, value_list

if __name__ == '__main__':
    tmp = ParseReaxys()
    tmp.forward()