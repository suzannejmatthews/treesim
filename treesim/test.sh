#150: .341 857
#576: .5178 .9256
#525  .9005 .9924

for i in 1 2 3 4 5
do

./treesim 567 33306 .5178 .5178 -o out$i.tre -p
done

ls out*.tre | xargs -Istr treezip str

