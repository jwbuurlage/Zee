#!/usr/bin/python3

import os

files = os.popen('ls */*.tar.gz').read()
files = files[:len(files) - 1]

files = files.split("\n")
for filename in files:
    os.system("tar -xzf " + filename)

for i in range(0, len(files)):
    files[i] = files[i][files[i].index('/') + 1:len(files[i]) - 7]

os.system("mkdir -p ../mtx")


for filename in files:
    os.system("mv " + filename + "/" + filename + ".mtx ../mtx")
    os.system("cp ../mtx/*.mtx ~/homework/thesis/zee/data/matrices")

print(files)
f = open("collection.txt", "w")
f.write(" ".join(files))
f.close()
