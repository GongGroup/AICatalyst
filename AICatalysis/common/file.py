import hashlib
import json
from collections import defaultdict
from json import JSONEncoder
from pathlib import Path

import yaml


def md5(name: str):
    return hashlib.md5(name.encode(encoding='utf-8')).hexdigest()


def ftemp(file: Path):
    """
    Add "_" in the prefix of file, e.g., chemical.json -> _chemical.json

    Args:
        file: file path

    Returns:
        temp: new file path
    """
    return file.parent / f"_{file.name}"


def fcopy(file: Path):
    """
    Add "_copy" in the suffix of file, e.g., chemical.json -> chemical_copy.json

    Args:
        file: file path

    Returns:
        temp: new file path
    """
    return file.parent / f"{file.stem}_copy{file.suffix}"


class HtmlIO(object):
    @staticmethod
    def read(file, encoding='utf-8'):
        with open(file, 'r', encoding=encoding) as f:
            content = f.read()
        return content


class JsonIO(object):

    @staticmethod
    def read(file, encoding='utf-8'):
        """
        Json file read func

        Args:
            file: json file
            encoding: file encoding method

        Returns:
            data: python object load from json file
        """
        with open(file, 'r', encoding=encoding) as f:
            data = json.load(f)
        return data

    @staticmethod
    def write(obj, file, indent=2, encoding='utf-8'):
        """
        Json file write func

        Args:
            obj: python object you want to write in json
            file: json file
            indent: indent depth, default: 2
            encoding: file encoding method, default: utf-8

        """
        with open(file, "w", encoding=encoding) as f:
            json.dump(obj, f, indent=indent)


class CatalystJSONEncoder(JSONEncoder):
    """
    Extend the default JSONEncoder, make it accept the <MCatalyst> instance
    """

    def default(self, o):
        from species import Reagent

        if isinstance(o, Reagent):
            return o.name
        else:
            return super(CatalystJSONEncoder, self).default(o)


class CatalystJsonIO(JsonIO):
    @staticmethod
    def write(obj, file, indent=2, encoding='utf-8'):
        with open(file, "w", encoding=encoding) as f:
            json.dump(obj, f, cls=CatalystJSONEncoder, indent=indent)


class YamlIO(object):
    yaml.warnings({'YAMLLoadWarning': False})

    @staticmethod
    def read(file, encoding='utf-8'):
        with open(file, 'r', encoding=encoding) as f:
            data = f.read()
        return yaml.safe_load(data)


class QMIO(object):

    @staticmethod
    def read1(file, encoding='utf-8'):
        with open(file, 'r', encoding=encoding) as f:
            data = f.read()
        parameters = {}
        for line in data.splitlines()[1:]:
            parameter = line.split()
            parameters[parameter[0]] = [float(parameter[1]), float(parameter[2]), float(parameter[3])]
        return parameters

    @staticmethod
    def read2(file, encoding='utf-8'):
        with open(file, 'r', encoding=encoding) as f:
            data = f.read()
        parameters = defaultdict(dict)
        for line in data.splitlines()[1:]:
            parameter = line.split()
            parameters[parameter[0]].update(
                {parameter[1]: [float(parameter[2]), float(parameter[3]), float(parameter[4])]})
        return parameters


if __name__ == '__main__':
    pass