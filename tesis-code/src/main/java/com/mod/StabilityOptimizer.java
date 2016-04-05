package com.mod;

import java.io.BufferedWriter;
import java.io.Console;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

import org.apache.commons.lang.ArrayUtils;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

//import scala.actors.threadpool.Arrays;
//import scala.reflect.generic.Trees.This;

public class StabilityOptimizer {
	
  public static void main(String[] args) throws IOException
  {
    boolean update;
    double modularity, resolution2, stability,meanVI;
    int  i, j, nClusters;
    int[] cluster;
    Network network;
    Random random;
    
    String inputFileName = args[0];
    String outputFileName = args[1];
    double resolution = 1;
    int nRandomStarts = Integer.parseInt(args[2]);
    int nIterations = Integer.parseInt(args[3]);
    long seed = Long.parseLong(args[4]);
    boolean printOutput = (Integer.parseInt(args[5]) > 0);
    boolean getPartition = (Integer.parseInt(args[6])>0);//if true means we only want to optimize for 1 markovTime
    double markovTime = (getPartition) ? Double.parseDouble(args[7]) : 1;
        
    BiMap<String,Integer> nodeBiMap = createMapping(inputFileName);
    System.gc();
    DoubleMatrix A = readFile(inputFileName, nodeBiMap);
    // DoubleMatrix A = readMappedFile(inputFileName);

        
    int n = nodeBiMap.size();
    // int n = 262111;
    double[] dataForOnes = new double[n];
    Arrays.fill(dataForOnes,1);
    DoubleMatrix ones = new DoubleMatrix(dataForOnes);
    DoubleMatrix d = A.mmul(ones); 
    double m = d.dot(ones)/2;
    DoubleMatrix pi = d.transpose().div(2*m);     	
    double[][] dataForPI = new double[n][n];    
    double[][] dataForD = new double[n][n];
    double[][] dataForDinv = new double[n][n];
    for (i=0; i<n; i++){
      dataForD[i][i] = d.get(i);
      dataForDinv[i][i]= 1/d.get(i);
      dataForPI[i][i]=pi.get(i); 
    }
    DoubleMatrix D = new DoubleMatrix(dataForD);
    DoubleMatrix Dinv = new DoubleMatrix(dataForDinv);
    DoubleMatrix L = D.sub(A);
    DoubleMatrix DinvL = Dinv.mmul(L);
    DoubleMatrix PI = new DoubleMatrix(dataForPI);
    DoubleMatrix piTpi = (pi.transpose()).mmul(pi);
        
    double[] times = new double[]{0.1,0.107189,0.110975,0.118953,0.123155,0.127505,0.132009,0.136672,0.146497,0.151672,0.162576,0.168318,0.174263,0.180419,0.186791,0.193389,0.20022,0.222195,0.238169,0.255291,0.264308,0.273644,0.293317,0.303677,0.314404,0.337006,0.361234,0.373994,0.387204,0.41504,0.4297,0.444878,0.460592,0.493705,0.511143,0.567243,0.587279,0.608022,0.629499,0.651734,0.674754,0.723263,0.74881,0.77526,0.802643,0.860346,0.890735,0.922198,0.954772,1.02341,1.05956,1.13573,1.17585,1.26038,1.3049,1.39871,1.44812,1.55223,1.60705,1.72259,1.84642,1.91164,1.97917,2.04907,2.12145,2.27397,2.35429,2.43744,2.52354,2.61268,2.70496,2.8005,2.89942,3.00184,3.21764,3.33129,3.44896,3.57079,3.69691,3.82749,3.96269,4.10266,4.24757,4.3976,4.55294,4.71375,4.88025,5.05263,5.2311,5.41587,5.60717,5.80523,6.01028,6.22257,6.44236,6.66992,6.90551,7.14943,7.40196,7.66341,7.9341,8.21434,8.50449,8.80488,9.11589,9.43788,9.77124,10.1164,10.4737,10.8437,11.2267,11.6232,12.0338,12.4588,12.8989,13.3545,13.8262,14.3146,14.8202,15.3437,15.8857,16.4468,17.0277,17.6291,18.2518,18.8965,19.564,20.255,20.9705,21.7112,22.4781,23.272,24.094,24.9451,25.8262,26.7384,27.6829,28.6607,29.673,30.7211,31.8063,32.9297,34.0929,35.2971,36.5438,37.8346,39.171,40.5546,41.9871,43.4701,45.0056,46.5953,48.2411,49.9451,51.7092,53.5357,55.4266,57.3844,59.4113,61.5099,63.6825,65.9319,68.2607,70.6718,73.1681,75.7525,78.4282,81.1984,84.0665,87.0359,90.1102,93.293,96.5883,100};
    // double[] times = new double[]{0.303677,0.511143,0.723263,0.922198,1.55223,3.00184,5.05263,7.14943,9.77124,12.0338,14.8202,18.2518,22.4781,30.7211,40.5546,51.7092,65.9319,81.1984,100};
    //double[] times = new double[]{0.314404,0.41504,0.511143,0.601734,0.802643,0.914772,1,2.04907,3.00184,10.8437,20.255,30.7211,40.5546,61.5099,100};
	
    ArrayList<Cluster> clusters = new ArrayList<Cluster>(); 
    ArrayList<Cluster> bestClusters = new ArrayList<Cluster>(); 
        
        
    if(getPartition){
      times = new double[]{markovTime};
    }
    
    for (double time: times){
      network = createNetwork(DinvL, PI, n, ones,time);
      if (printOutput)
        {
          System.out.format("Number of nodes: %d%n", network.getNumNodes());
          System.out.format("Number of edges: %d%n", network.getNumEdges()/2);
          System.out.println();
          System.out.println("Running Louvain algorithm");
          System.out.println();
        }
      
      resolution2 = resolution / network.getTotalEdgeWeight();
      cluster = null;
      nClusters = -1;
      random = new Random(seed);
      System.out.println("Starting random starts");
      for (i = 0; i < nRandomStarts; i++)
        {
          network.initSingleClusters();
          j = 0;
          update = true;
          do
            { 
              update = network.runLouvain(resolution2, random);  
              j++;
            }
          while ((j < nIterations) && update);
          modularity = network.calcModularity(resolution2);
          network.orderClustersByNumNodes();
          cluster = network.getClusters();
          nClusters = network.getNumClusters();
          System.out.println("Number of clusters found: " + nClusters);
          DoubleMatrix H = getPartitionMatrix(cluster,n, nClusters);
          stability = computeStability(H,PI, piTpi, DinvL,time);
          clusters.add(new Cluster (cluster, stability, nClusters,time));
          if (nClusters < 2){
            System.out.println("a.break!");
            break;
          }
        }
      System.out.println("Finish random starts. Finiding max stability cluster");
      
      //HACER QUE c SEA EL CLUSTER DE MAXIMA ESTABILIDAD
      Cluster best = clusters.get(0);
      for ( i = 1; i < clusters.size(); i++)
        {
          if (best.stability < clusters.get(i).stability) 
            best = clusters.get(i);
        }
      meanVI = computeMeanVI(clusters);
      best.setMeanVI(meanVI);
      bestClusters.add(best);
      clusters.clear();
      if (nClusters < 3){
        System.out.println("b.break!");
        break;
      }
    }
        
    if(getPartition){
      bestClusters.get(0).writeCluster(outputFileName,nodeBiMap);
    }else{
      writeResults(outputFileName, bestClusters);
    }
        	
  }
  
  private static BiMap<String,Integer> createMapping(String fileName) throws IOException{
    Scanner fileScanner = new Scanner(new FileReader(fileName));
    ArrayList<String> nodeList = new ArrayList<String>();
    fileScanner.nextLine(); //to supress headers
    while (fileScanner.hasNext())
      {
        Scanner lineScanner = new Scanner(fileScanner.nextLine());
        String node1 = lineScanner.next();
        String node2 = lineScanner.next();
        if (!node1.equals(node2))
          {
            nodeList.add(node1);
            nodeList.add(node2);
          }
      }
    fileScanner.close();
    System.gc();   
    Set<String> nodeSet = new HashSet<String>(nodeList);
    List<String> sortedNodes = new ArrayList<String>(nodeSet);
    Collections.sort(sortedNodes);
    BiMap<String, Integer> nodeBiMap = HashBiMap.create();
    int i;
    for (i = 0; i < sortedNodes.size(); i++ ){
      nodeBiMap.put(sortedNodes.get(i), i);
    }
    return nodeBiMap;    	 	
  }
  

  private static DoubleMatrix readMappedFile (String fileName) throws IOException{ 
    Scanner fileScanner, lineScanner;
    // double weight;
    int node1, node2;
    int n = 262111;
    
    // double[][] dataForA = new double[n][n];
    System.out.println(n);

    DoubleMatrix A = DoubleMatrix.zeros(n,n);
    
    fileScanner = new Scanner(new FileReader(fileName));
    while (fileScanner.hasNext())
      {
        lineScanner = new Scanner(fileScanner.nextLine());
        node1 = Integer.parseInt(lineScanner.next());
        node2 = Integer.parseInt(lineScanner.next());
        A.put(node1,node2, 1);
        A.put(node2,node1,1);
       
        // dataForA[node1][node2] = 1;
        // dataForA[node2][node1] = 1;
      } 
    // DoubleMatrix A = new DoubleMatrix(dataForA);
    return A;
  }


  private static DoubleMatrix readFile (String fileName, BiMap<String,Integer> nodeBiMap) throws IOException{	
    Scanner fileScanner, lineScanner;
    double weight;
    int node1, node2;
    int n = nodeBiMap.size();
    System.out.println("Node bimap size");
    System.out.println(n);
    
    double[][] dataForA = new double[n][n];
    
    fileScanner = new Scanner(new FileReader(fileName));
    fileScanner.nextLine();
    while (fileScanner.hasNext())
      {
        lineScanner = new Scanner(fileScanner.nextLine());
        node1 = nodeBiMap.get(lineScanner.next());
        node2 = nodeBiMap.get(lineScanner.next());
        weight = lineScanner.hasNextDouble() ? lineScanner.nextDouble() : 1;
        dataForA[node1][node2] = weight;
        dataForA[node2][node1] = weight;
      } 
    DoubleMatrix A = new DoubleMatrix(dataForA);
    return A;
  }
  
  private static Network createNetwork(DoubleMatrix DinvL, DoubleMatrix PI, int n, DoubleMatrix ones, double time) throws IOException{
    
    double weight;
    int node1, node2;
    DoubleMatrix exponent = DinvL.mul(-1*time);
    DoubleMatrix Pt = MatrixFunctions.expm(exponent); 
    DoubleMatrix At = (PI.mmul(Pt)); 	
    
    double[][] matrix = At.toArray2();
    ArrayList<Edge> edgeArrayList = new ArrayList<Edge>();
    
    for (int i =0; i<n;i++){
      for (int j=0; j<n; j++){
        if(i!=j){
          weight = matrix[i][j];
          if(!(weight < 0.00000001 && weight > -0.00000001)){
            edgeArrayList.add(new Edge(i,j,weight));
            edgeArrayList.add(new Edge(j,i,weight));
          }		
        }
      }
    }
    int numEdge = edgeArrayList.size();
    //System.out.println("numEdge:"+numEdge);
    
    Collections.sort(edgeArrayList);
    
    ArrayList<Edge> edgeArrayList2 = new ArrayList<Edge>();
    node1 = -1;
    node2 = -1;
    double edgeWeight = 0;
    for (int i = 0; i < numEdge; i++)
      {
        Edge edge = edgeArrayList.get(i);
        if ((edge.node1 == node2) && (edge.node2 == node1)||(edge.node1 == node1) && (edge.node2 == node2))
          edgeWeight += edge.weight;
        else
          {
            if (i > 0)
              edgeArrayList2.add(new Edge(node1, node2, edgeWeight));
            node1 = edge.node1;
            node2 = edge.node2;
            edgeWeight = edge.weight;
          }
      }
    edgeArrayList2.add(new Edge(node1, node2, edgeWeight));
    
    int numEdge2 = edgeArrayList2.size();
    int[][] edge2= new int[numEdge2][2];
    double[] edgeWeight2 = new double[numEdge2];
    for (int i=0; i<numEdge2; i++){
      Edge edge = edgeArrayList2.get(i);
      edge2[i][0] = edge.node1;
      edge2[i][1] = edge.node2;
      edgeWeight2[i] = edge.weight;
    }
    
    //System.out.println("edge2.size:"+edge2.length);
    //System.out.println("edgeWeight2.size" + edgeWeight2.length);
    double[] nodeWeight = At.mmul(ones).toArray();
    Network network = new Network(n, edge2, edgeWeight2, nodeWeight);
    
    return network;
  }
  
  private static void writeResults(String fileName, ArrayList<Cluster> bestClusters) throws IOException{
    BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(fileName));
    bufferedWriter.write("MarkovTime\tStability\tnumClusters\tmeanVI");
    bufferedWriter.newLine();
    for (int i = 0; i<bestClusters.size();i++){
      bufferedWriter.write(bestClusters.get(i).toString());
      bufferedWriter.newLine();
    }	
    bufferedWriter.close();    
  }
  
  /*private static void writeCluster(String fileName, Cluster cluster, BiMap<String,Integer> nodeBiMap) throws IOException
  {
    BufferedWriter bufferedWriter;
    int i;
    String word;    
    Integer index;
    
    bufferedWriter = new BufferedWriter(new FileWriter(fileName));

    int[] partition = cluster.partition;
    int numCommunities = cluster.numCommunities;

    //System.out.println(nodeBiMap);
    //System.out.println(Arrays.toString(partition));
    
    for (i = 0; i < partition.length; i++)
      {
        word = nodeBiMap.inverse().get(i);
        bufferedWriter.write(word + "\t"+Integer.toString(partition[i]));
        bufferedWriter.newLine();
      }
    bufferedWriter.close();        
  }*/
  
  private static double computeStability(DoubleMatrix H, DoubleMatrix PI, DoubleMatrix piTpi, DoubleMatrix DinvL, double time){
    DoubleMatrix exponent = DinvL.mul(-1*time);
    DoubleMatrix Pt = MatrixFunctions.expm(exponent);
    DoubleMatrix inside = (PI.mmul(Pt)).sub(piTpi);
    DoubleMatrix Rt = (H.transpose()).mmul(inside).mmul(H);
    int n = Rt.columns;
    double trace = 0;
    for (int i=0; i<n; i++){
      trace +=Rt.get(i,i);
    }
    return 0.5*trace;
  }

  private static DoubleMatrix getPartitionMatrix(int[] cluster, int n, int nClusters){
    DoubleMatrix H = DoubleMatrix.zeros(n,nClusters);
    for (int i=0; i< cluster.length; i++)
      H.put(i, cluster[i],1);    		
    return H;
  }
	
  private static double VI(Cluster c1, Cluster c2){
    double pk1, pk2, ent1=0, ent2=0, pkk, minfo=0, vi;
    int[] part1 = c1.partition;
    int[] part2 = c2.partition;
    int m1 = c1.numCommunities;
    int m2 = c2.numCommunities;
    int n = part1.length;
    DoubleMatrix H1 = getPartitionMatrix(part1, n, m1 );
    DoubleMatrix H2 = getPartitionMatrix(part2, n, m2);
    DoubleMatrix Pk1 = (H1.transpose()).mmul(H1);
    DoubleMatrix Pk2 = (H2.transpose()).mmul(H2); 	
    DoubleMatrix Pkk = (H1.transpose()).mmul(H2);
    for (int i =0; i< m1 ; i++){
      pk1 = Pk1.get(i,i)/n;
      if(Math.abs(pk1)>0.000001 || pk1 > 0){
        ent1 = ent1 + pk1*Math.log(pk1);
      }
      for ( int j = 0; j<m2 ; j++){
        pk2 = Pk2.get(j,j)/n;
        if ( i==0){
          if(Math.abs(pk1)>0.000001 || pk1 > 0){
            ent2 = ent2 + pk2*Math.log(pk2);
          }
        }
        pkk = Pkk.get(i,j)/n;
        if((Math.abs(pkk)>0.000001 || pkk > 0) && Math.abs(pk1*pk2)>0.000001){
          minfo = minfo + (pkk*Math.log(pkk/(pk1*pk2)));	
        }
      }
    		
    } 
    ent1 = -1*ent1;
    ent2 = -1*ent2;
    vi = ent1 + ent2 - 2*minfo;
    //System.out.println(" ;ent1: "+ent1+ " ;ent2: "+ent2+ ";minfo: "+minfo+" ;VI: "+vi);
    return vi;
  }
    
  private static double computeMeanVI(ArrayList<Cluster> clusters){
    double sumVI = 0;
    Cluster c1;
    int m = clusters.size();
    for (int i = 0; i< m; i++){
      c1 = clusters.get(0);
      clusters.remove(0);
      for (int j=0 ; j < m-(i+1); j++)
        sumVI = sumVI + VI(c1, clusters.get(j));		
    }
    return sumVI/(m*(m-1)/2);	
  }
	
}
