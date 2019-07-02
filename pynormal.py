import pynormal
import numpy as np

nodeNum = 249                 # 网格节点个数
triNum = 494                  # 三角形的个数
start = 0

nodes = np.zeros(nodeNum * 3, dtype='double')
trias = np.zeros(triNum * 3, dtype='int')

with open('points.txt') as f:
    i = 0
    for line in f:
        data = line.split(' ')
        nodes[i*3+0] = float(data[1])
        nodes[i*3+1] = float(data[2])
        nodes[i*3+2] = float(data[3])
        i += 1

with open('triangles.txt') as f:
    i = 0
    for line in f:
        data = line.split(' ')
        trias[i*3+0] = int(data[0])
        trias[i*3+1] = int(data[1])
        trias[i*3+2] = int(data[2])
        i += 1


result = pynormal.pynormal(nodes, trias)

for i in range(nodeNum):
    print(1-(nodes[i*3+0]*result[i*3+0]+nodes[i*3+1]*result[i*3+1]+(nodes[i*3+2]-0.005)*result[i*3+2])/0.002)
    