from bs4 import BeautifulSoup
import json

def parse_RX(reaction_type, type_index, RX_config):
    key_list = []
    value_list = []
    
    for event, elem in RX_config.items():
        key_list.append(event)
    
        elem_text = reaction_type.find_all(elem)
        elem_text_list = []
        for i in range(len(elem_text)):
            elem_text_list.append(elem_text[i].text)
        value_list.append(elem_text_list)
    
    list2dict(key_list, value_list, type_index)

def parse_RXD(reaction_type, type_index, RXD_config, RXD_index):
    cur_RXD = RXD_list[RXD_index]
    key_list = []
    value_list = []
        
    for event, elem in RXD_config.items():
        key_list.append(event)
        
        elem_value_list = []
        elem_list = cur_RXD.find_all(elem)
        for i in range(len(elem_list)):
            elem_value_list.append(elem_list[i].text)
        value_list.append(elem_value_list)
    list2dict(key_list, value_list, type_index)

def list2dict(key_list, value_list, type_index):
    dict_reaction = dict(zip(key_list, value_list))
    data[type_index].append(dict_reaction) 

if __name__ == '__main__':
    with open('RX_config.json', 'r') as f:
        RX_config = json.load(f)
    
    with open('RXD_config.json', 'r') as f:
        RXD_config = json.load(f)
    
    with open('reaxys_xml.xml', 'r') as f:
        xml_file = f.read()

    reaction_etree = BeautifulSoup(xml_file, 'xml')
    reaction_type_list = reaction_etree.find_all(name='reaction')
    reaction_type_num = len(reaction_type_list)
    data = [[] for _ in range(reaction_type_num)]

    for type_index in range(0, reaction_type_num):
        reaction_type = reaction_type_list[type_index]

        parse_RX(reaction_type, type_index, RX_config)
        
        RXD_list = reaction_type.find_all('RXD')
        for RXD_index in range(len(RXD_list)):
            parse_RXD(reaction_type, type_index, RXD_config, RXD_index)

    with open('../chemical/reaxys.json', 'w') as f:
        json.dump(data, f)
