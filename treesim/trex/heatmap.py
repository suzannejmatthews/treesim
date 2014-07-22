import matplotlib.pyplot as plt
import numpy as np
import sys

if __name__ == "__main__":
  if len(sys.argv)!= 3:
    print "usage: python", sys.argv[0], "<vals> <taxa>"
    print "where <vals> is the file that has all distances"
    print "and <taxa> is the number of taxa"
    sys.exit()

  column_labels = ['run1','run2','run3','run4','run5']
  row_labels = ['run1','run2','run3','run4','run5']
  data = np.random.rand(5,5)

  taxa = int(sys.argv[2])
  myF = open(sys.argv[1], "r")

  i,j = 0,0
  for line in myF:
    val = float(line.strip())/taxa
    data[i][j] = val
    j+=1
    if j == 5:
      j = 0
      i+=1 

  print data
  fig, ax = plt.subplots()
  heatmap = ax.pcolor(data, cmap=plt.cm.hot, vmin=0, vmax=1)

  # put the major ticks at the middle of each cell
  ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
  ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

  # want a more natural, table-like display
  ax.invert_yaxis()
  ax.xaxis.tick_top()

  ax.set_xticklabels(row_labels, minor=False)
  ax.set_yticklabels(column_labels, minor=False)
  plt.show()
