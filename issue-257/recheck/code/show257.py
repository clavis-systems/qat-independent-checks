import json
d=json.load(open('/s/issue_257.json')); b=d['body']
print(b[:6500])
