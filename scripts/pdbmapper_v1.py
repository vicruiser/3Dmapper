### PDBMapper (provisional name)
print("hello")


import re
with open("filename") as origin_file:
    for line in origin_file:
        line = re.findall(r'something', line)
        if line:
            line = line[0].split('"')[1]
        print ("line")