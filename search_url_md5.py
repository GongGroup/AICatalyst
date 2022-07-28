import json

with open("md5_name.json", "r", encoding="utf-8") as f:
    md5_dict = json.load(f)

print(md5_dict['8b97fdb93f66274f868e26ce4c5835c6'])