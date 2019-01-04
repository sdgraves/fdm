import subprocess
import os
import sys
import time
import mmap
import struct
import numpy as np
import matplotlib.pyplot as pyplot

fd = int(sys.argv[1])
num = int(sys.argv[2])
print("Starting python with fd ", fd)
print("mmapping length ", num)
mm = mmap.mmap(fd, num)
#size should be 24
header_format = "ddQ"
length = 1
off = 0
data = []
while not (length == 0):
    (T, t0, length ) = struct.unpack_from( header_format, mm, offset=off )
    off += 24
    print("Found: ", length, ' ', T, ' ', t0)
    d = struct.unpack_from(str(length) + 'd', mm, offset=off)
    data.append(d)
    off += length*8
    print("Next read at offset: ", off )
    
mm.close()

pyplot.plot(data[0])
#print(data[0])

os.close(fd)
print("DONE")
pyplot.show()
