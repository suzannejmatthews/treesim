# treesim
TreeSim: Accurate Simulation of Large Collections of Phylogenetic Trees

This is the project page for TreeSim, an efficient software package that enables the simulation of large collections of phylogenetic trees from existing consensus data. The software package contains three subdirectories. The `unresolverstr` and `unresolvermaj` directories contain my *fixed* versions of the original software provided by Sul and Williams. The `treesim` folder contains the treesim code that I developed. I am very grateful to Dr. Seung-Jin Sul and Dr. Tiffani Williams for making `unresolverstr` and `unresolvermaj` available. 

## Using TreeSim
TreeSim is very easy to use, and has many different options. You can generate unweighted and weighted trees, and also specify consensus trees or simply consensus rates. The algorithms you can use to generate phylogenetic tree collections include the combined consensus, the strict consensus, and majority consensus. TreeSim also allows you quickly generate random trees. The `-w` flag will allow you to simulate weighted trees.

###Simulating Collections from Consensus Trees
Use the `-m` and `-s` flags to specify the majority and strict consensus tree files. TreeSim will run the appropriate algorithm based on the combination of flags you provide.

To generate simulate a collection of 567 taxa and 3177 trees from a majority consensus tree specified in `run.0.maj`, run the following command:

`./treesim 567 3177 0 0 -m run.0.maj -o 567maj.0`

This generates 3,177 567-taxa trees using the majority simulation algorithm and stores the result in output file `567maj.0`.  

Similarily, the following command allows you to simulate a collection from a strict consensus tree:

`./treesim 567 3177 0 0 -s run.0.str -o 567str.0`

This will use the strict consensus simulation algorithm to produce the trees that will be stored in `567str.0`.

Running the combined consensus algorithm requires both a strict and majority consensus tree. Of all avaliable approaches, this algorithm produces tree collections that have a level of topological diversity that most closely resembles that found in real collections produced by phylogenetic search. The command:
`./treesim 567 3177 0 0 -s run.0.str -m run.0.maj -o 567comb.0`
will simulate a collection using the combined consensus algorithm.

###Simulating Collections from Consensus Rates
If the `-m` and `-s` flags are not used, TreeSim will ask that you provide the strict and/or majority consensus rates to simulate a collection. For example:

`./treesim 567 3177 .66 .66 -o 567str`

generates a tree collecton with strict consensus resolution of 66%. 

In contrast, the command:

`./treesim 567 3177 0 .9531 -o 567maj`

generates a tree collection with majority consensus resolution 95.31%.

Providing both a strict and majority consensus rates will simulate a collection using the combined consensus algorithm:

`./treesim 567 3177 .66 .9531 -o 567comb`

The trees produced by the above command will have strict resolution of 66% and a majority consensus resolution of 95.31%. 

###Simulating random trees
If you don't have any consensus trees or consensus rates, your only option is to use the random option. A collection 3,177 random 567-taxa trees can be generated with the following command: 

`./treesim 567 3177 0 0 -o 567rand`




## Available Consensus Trees
We make the strict and majority consensus trees generated from HashCS available in the `consensus` subdirectory. This should allow you to simulate collections similar to those described in our paper.


##Citation info
If you use TreeSim, we ask that you cite the following publication:

Matthews SJ. "Accurate Simulation of Large Collections of Phylogenetic Trees". In Proceedings of the 2015 IEEE International Conference on Bioinformatics and Biomedicine (BIBM'15). To appear.
