"""
Config reader.
Since the main script is split in two parts, a config makes it easy to \n
reuse the same settings.
"""

def readConfig(config):

    confDict = {}

    for line in config:
        if "=" in line:
            v = line.split("=")
            confDict[v[0]] = v[1].strip()

    return(confDict)

# if __name__ == '__main__':
#     with open("../config.txt") as cf:
#         config = readConfig(cf.readlines())
