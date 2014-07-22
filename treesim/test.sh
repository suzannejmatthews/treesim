for i in 1 2 3 4 5
do
./treesim 567 33306 0 .9256 -o out$i.tre
done

ls out*.tre | xargs -Istr treezip str

