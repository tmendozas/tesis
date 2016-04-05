package com.mod;

import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedWriter;

import com.google.common.collect.BiMap;

public class Cluster implements Comparable<Cluster>{
	
	public int[] partition;
	public double stability;
	public int numCommunities;
	public double resolution; //the resolution used to obtain this cluster (ex. markov time)
	public double meanVI; //not necessary, but used for best clusters 
	
	public Cluster(int[] partition, double stability, int numCommunities, double resolution) {
		this.partition = partition;
	    this.stability = stability;
	    this.numCommunities = numCommunities;
	    this.resolution = resolution;
	}
	
	public Cluster(int[] partition, double stability, int numCommunities, double resolution, double meanVI) {
		this(partition, stability, numCommunities, resolution);
		this.meanVI = meanVI;
	}

	public int compareTo(Cluster cluster) {
		int ret = 0;
		if (stability != cluster.stability){
			if (stability > cluster.stability)
				ret = 1;
			else
				ret = -1;
		}
		return ret;	
	}
	
	public String toString(){
		String str = resolution+ "\t" + stability + "\t" + numCommunities + "\t" + meanVI;
		return str;
	}

	public void setMeanVI(double meanVI){
		this.meanVI = meanVI;
	}

	public void writeCluster(String fileName, BiMap<String,Integer> nodeBiMap) throws IOException
  {
    BufferedWriter bufferedWriter;
    int i;
    String word;    
    Integer index;
    
    bufferedWriter = new BufferedWriter(new FileWriter(fileName));

    int[] partition = this.partition;
    int numCommunities = this.numCommunities;
    for (i = 0; i < partition.length; i++)
      {
        word = nodeBiMap.inverse().get(i);
        bufferedWriter.write(word + "\t"+Integer.toString(partition[i]));
        bufferedWriter.newLine();
      }
    bufferedWriter.close();        
  }

	public void writeMappedCluster(String fileName) throws IOException
  {
    BufferedWriter bufferedWriter;
    int i;
    String word;    
    Integer index;
    
    bufferedWriter = new BufferedWriter(new FileWriter(fileName));

    int[] partition = this.partition;
    int numCommunities = this.numCommunities;
    for (i = 0; i < partition.length; i++)
      {
        bufferedWriter.write(i + "\t"+Integer.toString(partition[i]));
        bufferedWriter.newLine();
      }
    bufferedWriter.close();        
  }


}
