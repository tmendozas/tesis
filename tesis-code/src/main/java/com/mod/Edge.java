package com.mod;

public class Edge implements Comparable<Edge> {
    public int node1;
    public int node2;
    public double weight;

    public Edge(int node1, int node2, double weight) {
       this.node1 = node1;
       this.node2 = node2;
       this.weight = weight;
    }

    public int compareTo(Edge edge) {
       return (node1 == edge.node1) ? (node2 - edge.node2) : (node1 - edge.node1);
    }
    
    public String toString(){
    	String str = "("+node1+","+node2+","+weight+")";
		return str;
    	
    }
    
}
