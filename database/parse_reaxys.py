from bs4 import BeautifulSoup
import json

data_type = []
data_react = []


def list2dict(func):
    def decorate(*args, **kwargs):
        if func.__name__ == 'parse_RX':
            key_list, value_list = func(*args, **kwargs)
            dict_type = dict(zip(key_list, value_list))
            data_type.append(dict_type)
        elif func.__name__ == 'parse_RXD':
            key_list, value_list = func(*args, **kwargs)
            dict_react = dict(zip(key_list, value_list))
            dict_react.update(data_type[-1])
            data_react.append(dict_react)

    return decorate


@list2dict
def parse_RX(reaction_type, RX_config):
    key_list = []
    value_list = []

    for event, elem in RX_config.items():
        key_list.append(event)

        elem_text = reaction_type.find_all(elem)
        elem_text_list = []
        for i in range(len(elem_text)):
            elem_text_list.append(elem_text[i].text)
        value_list.append(elem_text_list)
    return key_list, value_list


@list2dict
def parse_RXD(RXD_config, cur_RXD):
    key_list = []
    value_list = []

    for event, elem in RXD_config.items():
        key_list.append(event)

        elem_value_list = []
        elem_list = cur_RXD.find_all(elem)
        for i in range(len(elem_list)):
            elem_value_list.append(elem_list[i].text)
        value_list.append(elem_value_list)
    return key_list, value_list


if __name__ == '__main__':
    with open('RX_config.json', 'r') as f:
        RX_config = json.load(f)

    with open('RXD_config.json', 'r') as f:
        RXD_config = json.load(f)

    with open('reaxys_xml.xml', 'r') as f:
        xml_file = f.read()

    reaction_etree = BeautifulSoup(xml_file, 'xml')
    reaction_type_list = reaction_etree.find_all(name='reaction')

    for reaction_type in reaction_type_list:

        parse_RX(reaction_type, RX_config)

        RXD_list = reaction_type.find_all('RXD')
        for cur_RXD in RXD_list:
            parse_RXD(RXD_config, cur_RXD)

    with open('reaxys_json.json', 'w') as f:
        json.dump(data_react, f)
