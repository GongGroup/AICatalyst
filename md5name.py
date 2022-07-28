import hashlib
import json

md5_dict = {}
with open("datafile.1.csv", encoding='utf-8') as f:
    content = f.readlines()

for line in content:
    url, doi = line.strip().split(",")[-2:]
    md5_name = hashlib.md5(url.encode(encoding='utf-8')).hexdigest()
    md5_dict[md5_name] = url

with open("md5_name.json", "w", encoding="utf-8") as f:
    json.dump(md5_dict, f)
