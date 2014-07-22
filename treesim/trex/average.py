import sys

if __name__ == "__main__":
  if len(sys.argv)!= 2:
    print "usage:", sys.argv[0], "<file>"
    print "where <file> is the matrix you want to average"
    sys.exit()

  myF = open(sys.argv[1], "r")
  total = 0.0
  count = 0
  for line in myF:
    line = line.strip().split()
    line = [int(i) for i in line]
    count += len(line)
    total += sum(line)

  myF.close()
  print total/count
