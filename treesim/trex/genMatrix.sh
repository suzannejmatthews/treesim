#!/bin/bash
count=1
path=~/Downloads/hashrf-pq-6.0
data=~/Documents/Research/TAMU/Combined/treefiles/567-taxa
sizes=(3177 3178 3178 3179 2985 2985 2986 2987 2162 2162 2163 2164)
for i in 0 1 2 3 4 5 6 7 8 9 10 11
do
  for j in 0 1 2 3 4 5 6 7 8 9 10 11
do
  echo "$path/hashrf $data/run.$i $data/run.$j ${sizes[i]} ${sizes[j]} -o o$count.txt"
  $path/hashrf $data/run.$i $data/run.$j ${sizes[i]} ${sizes[j]} -o o$count.txt
  python average.py o$count.txt >> final.matrix
  let count=count+1
done
done
