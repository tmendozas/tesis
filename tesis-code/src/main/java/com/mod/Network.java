package com.mod;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;

public class Network implements Serializable {

  private static final long serialVersionUID = 1;

  private int numNodes;
  private int[] firstNeighborIndex;
  private int[] neighbor;
  private double[] edgeWeight;
  private double totalEdgeWeightSelfLinks;
  private double[] nodeWeight;
  private int numClusters;
  private int[] cluster;

  private double[] clusterWeight;
  private int[] numNodesPerCluster;
  private int[][] nodePerCluster;
  private boolean statsAvailable;


  public Network(int numNodes, int[][] edge, double[] edgeWeight, double[] nodeWeight) {
    this(numNodes, edge, edgeWeight, nodeWeight, null);
  }

  public Network(int numNodes, int[][] edge, double[] edgeWeight, double[] nodeWeight, int[] cluster) {
    double[] edgeWeight2;
    int i, j, nEdges, nEdgesWithoutSelfLinks;
    int[] neighbor;

    this.numNodes = numNodes;

    nEdges = edge.length;
    firstNeighborIndex = new int[numNodes + 1];
    if (edgeWeight == null) 
    {
      edgeWeight = new double[nEdges];
      for (i = 0; i < nEdges; i++)
        edgeWeight[i] = 1;
    }

    totalEdgeWeightSelfLinks = 0;

    neighbor = new int[nEdges];
    edgeWeight2 = new double[nEdges];
    i = 1;
    nEdgesWithoutSelfLinks = 0;
    for (j = 0; j < nEdges; j++)
      if (edge[j][0] == edge[j][1])
        totalEdgeWeightSelfLinks += edgeWeight[j];
      else
      {
        if (edge[j][0] >= i)
          for (; i <= edge[j][0]; i++)
            firstNeighborIndex[i] = nEdgesWithoutSelfLinks;
        neighbor[nEdgesWithoutSelfLinks] = edge[j][1];
        edgeWeight2[nEdgesWithoutSelfLinks] = edgeWeight[j];
        nEdgesWithoutSelfLinks++;
      }
    for (; i <= numNodes; i++)
      firstNeighborIndex[i] = nEdgesWithoutSelfLinks;

    this.neighbor = new int[nEdgesWithoutSelfLinks];
    System.arraycopy(neighbor, 0, this.neighbor, 0, nEdgesWithoutSelfLinks);
    this.edgeWeight = new double[nEdgesWithoutSelfLinks];
    System.arraycopy(edgeWeight2, 0, this.edgeWeight, 0, nEdgesWithoutSelfLinks);

    if (nodeWeight == null)
    {
      this.nodeWeight = new double[numNodes];
      for (i = 0; i < numNodes; i++)
        this.nodeWeight[i] = 1;
    }
    else
      this.nodeWeight = nodeWeight;

    setClusters(cluster);
  }

  public int getNumNodes() { return numNodes; }

  public int getNumEdges() { return neighbor.length; }

  public int getNumClusters() { return numClusters; }

  public int[] getClusters() { return cluster; }

  public double getTotalEdgeWeight()
  {
    double totalEdgeWeight;
    int i;

    totalEdgeWeight = totalEdgeWeightSelfLinks;
    for (i = 0; i < neighbor.length; i++)
      totalEdgeWeight += edgeWeight[i];

    return totalEdgeWeight;
  }

  public void setClusters(int[] cluster)
  {
    int i, j;

    if (cluster == null)
      numClusters = 0;
    else
    {
      i = 0;
      for (j = 0; j < numNodes; j++)
        if (cluster[j] > i)
          i = cluster[j];
      numClusters = i + 1;
    }
    this.cluster = cluster;

    deleteStats();
  }

  public void initSingleClusters()
  {
    int i;

    numClusters = numNodes;
    cluster = new int[numNodes];
    for (i = 0; i < numNodes; i++)
      cluster[i] = i;

    deleteStats();
  }


  public void mergeClusters(int[] newCluster)
  {
    int i, j, k;

    if (cluster == null)
      return;

    i = 0;
    for (j = 0; j < numNodes; j++)
    {
      k = newCluster[cluster[j]];
      if (k > i)
        i = k;
      cluster[j] = k;
    }
    numClusters = i + 1;

    deleteStats();
  }

  public Network getReducedNetwork()
  {
    double[] reducedNetworkEdgeWeight1, reducedNetworkEdgeWeight2;
    int i, j, k, l, m, reducedNetworkNEdges1, reducedNetworkNEdges2;
    int[] reducedNetworkNeighbor1, reducedNetworkNeighbor2;
    Network reducedNetwork;

    if (cluster == null)
      return null;

    if (!statsAvailable)
      computeStats();

    reducedNetwork = new Network();

    reducedNetwork.numNodes = numClusters;
    reducedNetwork.firstNeighborIndex = new int[numClusters + 1];
    reducedNetwork.totalEdgeWeightSelfLinks = totalEdgeWeightSelfLinks;
    reducedNetwork.nodeWeight = new double[numClusters];

    reducedNetworkNeighbor1 = new int[neighbor.length];
    reducedNetworkEdgeWeight1 = new double[edgeWeight.length];

    reducedNetworkNeighbor2 = new int[numClusters - 1];
    reducedNetworkEdgeWeight2 = new double[numClusters];

    reducedNetworkNEdges1 = 0;
    for (i = 0; i < numClusters; i++)
    {
      reducedNetworkNEdges2 = 0;
      for (j = 0; j < nodePerCluster[i].length; j++) {
        k = nodePerCluster[i][j];
        for (l = firstNeighborIndex[k]; l < firstNeighborIndex[k + 1]; l++) {
          m = cluster[neighbor[l]];
          if (m != i) {
            if (reducedNetworkEdgeWeight2[m] == 0) {
              reducedNetworkNeighbor2[reducedNetworkNEdges2] = m;
              reducedNetworkNEdges2++; }
            reducedNetworkEdgeWeight2[m] += edgeWeight[l]; 
          } else
            reducedNetwork.totalEdgeWeightSelfLinks += edgeWeight[l];
        }
        reducedNetwork.nodeWeight[i] += nodeWeight[k];
      }

      for (j = 0; j < reducedNetworkNEdges2; j++)
      {
        reducedNetworkNeighbor1[reducedNetworkNEdges1 + j] = reducedNetworkNeighbor2[j];
        reducedNetworkEdgeWeight1[reducedNetworkNEdges1 + j] = reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[j]];
        reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[j]] = 0;
      }
      reducedNetworkNEdges1 += reducedNetworkNEdges2;
      reducedNetwork.firstNeighborIndex[i + 1] = reducedNetworkNEdges1;
    }

    reducedNetwork.neighbor = new int[reducedNetworkNEdges1];
    reducedNetwork.edgeWeight = new double[reducedNetworkNEdges1];
    System.arraycopy(reducedNetworkNeighbor1, 0, reducedNetwork.neighbor, 0, reducedNetworkNEdges1);
    System.arraycopy(reducedNetworkEdgeWeight1, 0, reducedNetwork.edgeWeight, 0, reducedNetworkNEdges1);
    return reducedNetwork;
  }


  public double calcModularity(double resolution)
  {
    double qualityFunction, totalEdgeWeight;
    int i, j, k;

    if (cluster == null)
      return Double.NaN;

    if (!statsAvailable)
      computeStats();

    qualityFunction = totalEdgeWeightSelfLinks;
    totalEdgeWeight = totalEdgeWeightSelfLinks;
    for (i = 0; i < numNodes; i++)
    {
      j = cluster[i];
      for (k = firstNeighborIndex[i]; k < firstNeighborIndex[i + 1]; k++)
      {
        if (cluster[neighbor[k]] == j)
          qualityFunction += edgeWeight[k];
        totalEdgeWeight += edgeWeight[k];
      }
    }
    for (i = 0; i < numClusters; i++)
      qualityFunction -= clusterWeight[i] * clusterWeight[i] * resolution;
    qualityFunction /= totalEdgeWeight;
    return qualityFunction;
  }


  public boolean runLocalMoves(double resolution, Random random)
  {
    boolean update;
    double maxQualityFunction, qualityFunction;
    double[] clusterWeight, edgeWeightPerCluster;
    int bestCluster, i, j, k, l, nNeighboringClusters, nStableNodes, nUnusedClusters;
    int[] neighboringCluster, newCluster, numNodesPerCluster, nodeOrder, unusedCluster;

    if ((cluster == null) || (numNodes == 1))
      return false;

    update = false;

    clusterWeight = new double[numNodes];
    numNodesPerCluster = new int[numNodes];
    for (i = 0; i < numNodes; i++)
    {
      clusterWeight[cluster[i]] += nodeWeight[i];
      numNodesPerCluster[cluster[i]]++;
    }

    nUnusedClusters = 0;
    unusedCluster = new int[numNodes];
    for (i = 0; i < numNodes; i++)
      if (numNodesPerCluster[i] == 0)
      {
        unusedCluster[nUnusedClusters] = i;
        nUnusedClusters++;
      }

    nodeOrder = new int[numNodes];
    for (i = 0; i < numNodes; i++)
      nodeOrder[i] = i;
    for (i = 0; i < numNodes; i++)
    {
      j = random.nextInt(numNodes);
      k = nodeOrder[i];
      nodeOrder[i] = nodeOrder[j];
      nodeOrder[j] = k;
    }

    edgeWeightPerCluster = new double[numNodes];
    neighboringCluster = new int[numNodes - 1];

    nStableNodes = 0;
    i = 0;
    do
    {
      j = nodeOrder[i];

      nNeighboringClusters = 0;
      for (k = firstNeighborIndex[j]; k < firstNeighborIndex[j + 1]; k++)
      {
        l = cluster[neighbor[k]];
        if (edgeWeightPerCluster[l] == 0)
        {
          neighboringCluster[nNeighboringClusters] = l;
          nNeighboringClusters++;
        }
        edgeWeightPerCluster[l] += edgeWeight[k];
      }

      clusterWeight[cluster[j]] -= nodeWeight[j];
      numNodesPerCluster[cluster[j]]--;
      if (numNodesPerCluster[cluster[j]] == 0)
      {
        unusedCluster[nUnusedClusters] = cluster[j];
        nUnusedClusters++;
      }

      bestCluster = -1;
      maxQualityFunction = 0;
      for (k = 0; k < nNeighboringClusters; k++)
      {
        l = neighboringCluster[k];
        qualityFunction = edgeWeightPerCluster[l] - nodeWeight[j] * clusterWeight[l] * resolution;
        if ((qualityFunction > maxQualityFunction) || ((qualityFunction == maxQualityFunction) && (l < bestCluster)))
        {
          bestCluster = l;
          maxQualityFunction = qualityFunction;
        }
        edgeWeightPerCluster[l] = 0;
      }
      if (maxQualityFunction == 0)
      {
        bestCluster = unusedCluster[nUnusedClusters - 1];
        nUnusedClusters--;
      }

      clusterWeight[bestCluster] += nodeWeight[j];
      numNodesPerCluster[bestCluster]++;
      if (bestCluster == cluster[j])
        nStableNodes++;
      else
      {
        cluster[j] = bestCluster;
        nStableNodes = 1;
        update = true;
      }

      i = (i < numNodes - 1) ? (i + 1) : 0;
    }
    while (nStableNodes < numNodes);

    newCluster = new int[numNodes];
    numClusters = 0;
    for (i = 0; i < numNodes; i++)
      if (numNodesPerCluster[i] > 0)
      {
        newCluster[i] = numClusters;
        numClusters++;
      }
    for (i = 0; i < numNodes; i++)
      cluster[i] = newCluster[cluster[i]];

    deleteStats();

    return update;
  }


  public boolean runLouvain(double resolution, Random random)
  {
      boolean update, update2;
      Network reducedNetwork;

      if ((cluster == null) || (numNodes == 1))
        return false;

      update = runLocalMoves(resolution, random);

      if (numClusters < numNodes)
      {
        reducedNetwork = getReducedNetwork();
        reducedNetwork.initSingleClusters();

        update2 = reducedNetwork.runLouvain(resolution, random);

        if (update2)
        {
          update = true;
          mergeClusters(reducedNetwork.getClusters());
        }
      }

      deleteStats();
      return update;
  }



  private Network()
  {
  }

  public void orderClustersByNumNodes()
  {
    class ClusterSize implements Comparable<ClusterSize>
    {
      public int cluster;
      public double size;

      public ClusterSize(int cluster, double size)
      {
        this.cluster = cluster;
        this.size = size;
      }

      public int compareTo(ClusterSize cluster)
      {
        return (cluster.size > size) ? 1 : ((cluster.size < size) ? -1 : 0);
      }
    }

    ClusterSize[] clusterSize;
    int i;
    int[] newCluster;

    if (cluster == null)
      return;
    if (!statsAvailable)
      computeStats();

    clusterSize = new ClusterSize[numClusters];
    for (i = 0; i < numClusters; i++)
      clusterSize[i] = new ClusterSize(i, numNodesPerCluster[i]);

    Arrays.sort(clusterSize);
    newCluster = new int[numClusters];
    i = 0;
    do
    {
      newCluster[clusterSize[i].cluster] = i;
      i++;
    }
    while ((i < numClusters) && (clusterSize[i].size > 0));
    numClusters = i;
    for (i = 0; i < numNodes; i++)
      cluster[i] = newCluster[cluster[i]];

    deleteStats();
  }

  

  private void computeStats()
  {
    int i, j;

    clusterWeight = new double[numClusters];
    numNodesPerCluster = new int[numClusters];
    nodePerCluster = new int[numClusters][];

    for (i = 0; i < numNodes; i++)
    {
      clusterWeight[cluster[i]] += nodeWeight[i];
      numNodesPerCluster[cluster[i]]++;
    }

    for (i = 0; i < numClusters; i++)
    {
      nodePerCluster[i] = new int[numNodesPerCluster[i]];
      numNodesPerCluster[i] = 0;
    }

    for (i = 0; i < numNodes; i++)
    {
      j = cluster[i];
      nodePerCluster[j][numNodesPerCluster[j]] = i;
      numNodesPerCluster[j]++;
    }

    statsAvailable = true;
  }

  private void deleteStats()
  {
    clusterWeight = null;
    numNodesPerCluster = null;
    nodePerCluster = null;
    statsAvailable = false;
  }
}
