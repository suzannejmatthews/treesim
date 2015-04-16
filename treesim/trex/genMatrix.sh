#!/bin/bash
count=1
path=~/Downloads/hashrf-pq-6.0
data=~/Documents/treefiles/525
#sizes=(3177 3178 3178 3179 2985 2985 2986 2987 2162 2162 2163 2164)
sizes=(30000 30000 30000 30000 30000)
for i in 0 1 2 3 4
do
  for j in 0 1 2 3 4
do
  echo "$path/hashrf $data/run.$i $data/run.$j ${sizes[i]} ${sizes[j]} -o o$count.txt"
  $path/hashrf $data/run.$i $data/run.$j ${sizes[i]} ${sizes[j]} -o o$count.txt
  python average.py o$count.txt >> final.matrix
  let count=count+1
done
done
