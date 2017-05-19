This thesis has two main objectives:

1. Explore communities at different times (times are hardcoded in StabilityOptimizer.java)

2. Get the partition in communities at a given time.

Usually, one would run the exploratory component to choose relevant times.
Then get the partition for one or several of those times.

Project should be compiled with maven, and .jar added to CLASSPATH:
export CLASSPATH=export CLASSPATH=target/tesis-0.1-jar-with-dependencies.jar

An example to run the exploratory component would be as follows:
>> java -Xmx14g com.mod.StabilityOptimizer ~/network.tsv ~/results/results.txt 5 5 1 1 0
Where the arguments are:
1.- Network adjacency file separated by tabs
2.- File where output should be written to
3.- Number of Random Starts
4.- Number of Iterations
5.- seed
6.- Print output? 1 = yes, 0 = no
7.- 0 = run exploratory component

To get the partition for a desired time, the arguments are similar.
Note that at the end, instead of a 0 we write the desired time.
>> java -Xmx14g com.mod.StabilityOptimizer ~/network.tsv ~/results/results.txt 5 5 1 1 1.84642 


