import org.apache.spark.mllib.linalg.distributed.{CoordinateMatrix, MatrixEntry,RowMatrix}
import org.apache.spark.rdd.RDD

// import breeze.linalg._
// import breeze.numerics._


// def main(args: Array[String]){

//   val sc = new SparkContext("local", "Amazon", System.getenv(
//       "SPARK_HOME"), System.getenv("SPARK_CLASSPATH").split(":"))

  val textFile = sc.textFile("/Users/tania/Downloads/amazon0302-nohead.txt")

  textFile.take(2)

  val graph = textFile.map(_.split("\t")).map(x => (x(0).toLong,x(1).toLong))

  graph.take(2)
  graph.count

  // val m1 = DenseMatrix((1.0,2.0), (3.0,4.0))

  //ONLY USE PRODUCTS FROM 1 TO 50,000
  // val graphjava = graph.filter( x => x._1 < 5000 && x._2 < 5000)
  val graphjava = graph

  val preundirect = graphjava.map(x => if (x._2 < x._1) (x._2,x._1) else x)


  val undirected = preundirect.map(x => (x,1)).reduceByKey(_+_)

  val source = undirected.map(x => (x._1._1,1))
  val target = undirected.map(x => (x._1._2,1))
  val degree = source.union(target).reduceByKey(_+_)

  val degMax = degree.map(x =>  x._2).reduce((a,b) => Math.max(a,b))
  val degMin = degree.map(x =>  x._2).reduce((a,b) => Math.min(a,b))
  val degMean = degree.map(x => x._2.toDouble).mean()
  val (degb, degc) = degree.map(x => x._2.toDouble).histogram(10)

  val degone = degree.filter(x => x._2==1)

  val degnotone = degree.filter(x => x._2>1)
  val reduced = undirected.map(x => x._1).join(degnotone).map(x => (x._2._1,x._1)).join(degnotone).map(x => (x._2._1,x._1))

  val redsource = reduced.map(x => (x._1,1))
  val redtarget = reduced.map(x => (x._2,1))
  val reddegree = redsource.union(redtarget).reduceByKey(_+_)
  val reddegone = reddegree.filter(x => x._2==1)

  reduced.map(item => "%s\t%s".format(item._1,item._2)).saveAsTextFile("/Users/tania/Downloads/amazon-0302-reduced")


/*Get a resume of all categories*/
val metaFile = sc.textFile("/Users/tania/Downloads/amazon-meta-head.txt")

val category = metaFile.filter(
  x => x.startsWith("   |")
  ).flatMap(
  _.split("""\|""")
  ).filter(
  x => !x.startsWith("  ")).map(x => (x,1)).reduceByKey(_+_).filter(x => x._2>1)

val categoryWithLevel = metaFile.filter(
  x => x.startsWith("   |")
  ).map(
  _.split("""\|""")
  ).flatMap(
  x => for (i <- 0 until x.length) yield ((x(i),i),1)
  ).filter( x => !x._1._1.startsWith(" ")
  ).reduceByKey(_+_).map(x => (x._2,x._1)).sortByKey(ascending=false)

categoryWithLevel.map(item => "%s\t%s\t%s".format(item._2._1,item._2._2,item._1)).saveAsTextFile("/Users/tania/Downloads/amazon-meta-categories")

/*Count how many items there are in a community*/
val graphFile = sc.textFile("/Users/tania/Downloads/amazon0302-reduced.tsv")
val graph = graphFile.map(_.split("\t"))
val src = graph.map(x => (x(0),1))
val target = graph.map(x => (x(1),1))
val degree = src.union(target).reduceByKey(_+_)

val comFile = sc.textFile("/Users/tania/tesis-itam/tesis/amazon-results/node-com4c.tsv")
val coms = comFile.map(_.split("\t")).map(x => (x(1),1))
val comCount = coms.reduceByKey(_+_)

val nodeInfo = degree.join(coms).map({case(id,(degree,com)) => 
  (degree,(id,com))}).sortByKey(ascending=false).map(x => 
    "%s\t%s\t%s".format(x._2._2,x._2._1,x._1)).saveAsTextFile("/Users/tania/Desktop/amazonData/nodecom-info")

comCount.collect
//   val symmetric = undirected.map(x => (x._1._1,x._1._2)).union(undirected.map(x => (x._1._2,x._1._1)))

//   val builder = new CSCMatrix.Builder[Double](rows=262111, cols=262111)
//   val collected = symmetric.collect

//   for(i <- 0 until collected.length){
//     val tuple = collected(i)
//     builder.add(tuple._1.toInt,tuple._2.toInt, 1.0)
// }
//   // etc.
//   val A = builder.result()

//   // val entries: RDD[MatrixEntry] = symmetric.map(x => MatrixEntry(x._1,x._2,1))
//   // val mat: CoordinateMatrix = new CoordinateMatrix(entries)
//   // val A = mat.toRowMatrix()

//   val ones = DenseVector.ones[Double](262111)

//   val d = A * ones
//   val m = d dot ones

//   val pi = d.t :/ (2*m)
//   val D = diag(d)
//   val Dinv = diag( ones :/ d)
//   val PI = diag(pi)
//   val L = D - A 
//   val DinvL = Dinv * L
//   val piTpi = pi.t * pi
//   val times = new double[]{0.303677,0.41504,0.511143,0.612377,0.723263,0.77526,0.860346,0.922198,1.02341,1.55223,2.04907,3.00184,4.10266,5.05263,7.14943,9.11589,9.43788,9.77124,12.0338,14.3146,14.8202,18.2518,22.4781,30.7211,40.5546,51.7092,65.9319,81.1984,100};



//   undirected.count

//   undirected.map(item => "%s\t%s".format(item._1, item._2)).saveAsTextFile("/Users/tania/Downloads/amazon0302-undirected")

// }