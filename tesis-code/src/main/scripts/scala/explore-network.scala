import org.apache.spark.mllib.linalg.distributed.{CoordinateMatrix, MatrixEntry,RowMatrix}
import org.apache.spark.rdd.RDD

import breeze.linalg._
import breeze.numerics._

import org.jblas.DoubleMatrix


// Delete all reviews.. 4 spaces followed by numbers
// perl -ni.bak -e'print unless m/^\s*\d+/' amazon-meta-head.txt
//perl -ni.bak -e'print unless m/^\s*similar/' amazon-meta-head.txt
// In vi remove single quotes...
//perl -ni.bak -e'print unless m/^\s*ASIN/' amazon-meta-head.txt
//perl -ni.bak -e 's/\n\s\sdiscontinued product/ discontinued/g' amazon-meta-head.txt 

val sc = new SparkContext("local", "Amazon", System.getenv(
      "SPARK_HOME"), System.getenv("SPARK_CLASSPATH").split(":"))

  val textFile = sc.textFile("/Users/tania/Downloads/amazon0302-nohead.txt")

  textFile.take(2)

  val graph = textFile.map(_.split("\t")).map(x => (x(0).toLong,x(1).toLong))

  graph.take(2)
  graph.count

  // val m1 = DenseMatrix((1.0,2.0), (3.0,4.0))


  val preundirect = graph.map(x => if (x._2 < x._1) (x._2,x._1) else x)


  val undirected = preundirect.map(x => (x,1)).reduceByKey(_+_)

  val symmetric_matrix_data = undirected.map(x => (x._1._1,x._1._2)).union(undirected.map(x => (x._1._2,x._1._1)))

  val source = undirected.map(x => (x._1,1))
  val target = undirected.map(x => (x._2,1))

  val nodes = source.union(target)

  val node_degree = nodes.reduceByKey(_+_)

  val degMax = node_degree.map(x => x._2).reduce((a,b) => Math.max(a,b))
    val degMin = node_degree.map(x => x._2).reduce((a,b) => Math.min(a,b))
    val degMean = node_degree.map(x => x._2.toDouble).mean()
    val (degb, degc) = node_degree.map(x => x._2.toDouble).histogram(10)


////////////////////////////
//Find nodes with highest degree by community

  val textFile = sc.textFile("/Users/tania/tesis-itam/tesis/amazon-results/amazon0302red-edges.tsv")

  val graph = textFile.map(_.split("\t")).map(x => (x(0),x(1)))
  val source = graph.map(x => x._1)
  val target = graph.map(x => x._2)
  val degree = source.union(target).map(x => (x,1)).reduceByKey(_+_)

  val coms = sc.textFile("/Users/tania/Desktop/amazonnode-community4c.tsv"
      ).map(_.split("\t")).map(x => (x(0),x(1)))
  
  val info = degree.join(coms).map(x => (x._2._1,(x._1,x._2._2))) //(deg,(node,com))

  val com0 = info.filter(x => x._2._2=="0").sortByKey(false)
  val com1 = info.filter(x => x._2._2=="1").sortByKey(false)
  val com2 = info.filter(x => x._2._2=="2").sortByKey(false)
  val com3 = info.filter(x => x._2._2=="3").sortByKey(false)

