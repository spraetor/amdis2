import os, sys
import fnmatch
import re

root = 'src'
pattern = re.compile('[#]include [<"]([a-zA-Z0-9/_]+.h)[>"]')

for d in os.listdir(root):
    dir = os.path.join(root, d)
    if os.path.isdir(dir):
        print "> ", dir
        for f in os.listdir(dir):
            filename = os.path.join(dir, f)
            if os.path.isfile(filename):
                # print f
                # with open(filename, "r+") as fp:
                #     lines = fp.readlines()
                #     fp.seek(0)
                #     fp.truncate()
                #     for line in lines:
                #         m = pattern.match(line)
                #         if m:
                #             include_f = m.group(1)
                #             if os.path.exists(os.path.join(root, include_f)):
                #                 line = line.replace(".h", ".hpp")
                #         fp.write(line)
    
                if fnmatch.fnmatch(filename, '*.cc'):
                    os.rename(filename, filename.replace(".cc", ".cpp"))
                elif fnmatch.fnmatch(filename, '*.h'):
                    os.rename(filename, filename.replace(".h", ".hpp"))