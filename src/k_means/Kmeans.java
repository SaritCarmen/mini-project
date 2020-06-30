package k_means;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import primitives.Point3D;

/**
 * @author rinat
 *
 */

public class Kmeans {

    //Number of Clusters. This metric should be related to the number of points
    private int NUM_CLUSTERS = 3;    
    private int MAX_ITERATIONS = 100;
    //private static final int MIN_COORDINATE = 0;
   // private static final int MAX_COORDINATE = 10;

    private List<Point> points;
    private List<Cluster> clusters;

    public Kmeans() {
        this.points = new ArrayList<Point>();
        this.clusters = new ArrayList<Cluster>();        
    }
   /* public KMeans(List<Point> points) {
        this.points = points;
        this.clusters = new ArrayList<Point>();        
    }*/

    
    /**
	 * @return the clusters
	 */
	public List<Cluster> getClusters() {
		return clusters;
	}

	//Initializes the process
    public void init(List<Point> points) {
        //add Points
        this.points = points;

        //Create Clusters
        //Set Random Centroids
        for (int i = 0; i < NUM_CLUSTERS; i++) {
            Cluster cluster = new Cluster(i);
            Point centroid = Point.createRandomPoint(-100,100);
            cluster.setCentroid(centroid);
            this.clusters.add(cluster);
        }

        //Print Initial state
        //plotClusters();
    }

    private void plotClusters() {
        for (Cluster c :clusters) {
            c.plotCluster();
        }
    }

    //The process to calculate the K Means, with iterating method.
    public void calculate() {
        boolean finish = false;
        int iteration = 0;

        // Add in new data, one at a time, recalculating centroids with each new one. 
        while(!finish && iteration<MAX_ITERATIONS) {
            //Clear cluster state
            clearClusters();

            List<Point> lastCentroids = getCentroids();

            //Assign points to the closer cluster
            assignCluster();

            //Calculate new centroids.
            calculateCentroids();

            iteration++;

            List<Point> currentCentroids = getCentroids();

            //Calculates total distance between new and old Centroids
            double distance = 0;
            for(int i = 0; i < lastCentroids.size(); i++) {
                distance += Point.distance(lastCentroids.get(i),currentCentroids.get(i));
            }
            /*System.out.println("#################");
            System.out.println("Iteration: " + iteration);
            System.out.println("Centroid distances: " + distance);
            plotClusters();*/

            if(distance == 0) {
                finish = true;
            }
        }
    }

    private void clearClusters() {
        for(Cluster cluster : clusters) {
            cluster.clear();
        }
    }

    private List<Point> getCentroids() {
        List<Point> centroids = new ArrayList<Point>(NUM_CLUSTERS);
        for(Cluster cluster : clusters) {
            Point aux = cluster.getCentroid();
            Point point = new Point(aux.getGeometry());
            centroids.add(point);
        }
        return centroids;
    }

    private void assignCluster() {
        double max = Double.MAX_VALUE;
        double min = max; 
        int cluster = 0;                 
        double distance = 0.0; 

        for(Point point : points) {
            min = max;
            for(int i = 0; i < NUM_CLUSTERS; i++) {
                Cluster c = clusters.get(i);
                distance = Point.distance(point, c.getCentroid());
                if(distance < min){
                    min = distance;
                    cluster = i;
                }
            }
            point.setCluster_number(cluster);
            clusters.get(cluster).addPoint(point);
        }
    }

    private void calculateCentroids() {
        for(Cluster cluster : clusters) {
            double sumX = 0;
            double sumY = 0;
            double sumZ = 0;
            List<Point> listPoints = cluster.getPoints();
            int n_points = listPoints.size();

            for(Point point : listPoints) {
                sumX += point.getPositionPoint().get_x().get();
                sumY += point.getPositionPoint().get_y().get();
                sumZ += point.getPositionPoint().get_z().get();}

            Point centroid = cluster.getCentroid();
            if(n_points > 0) {
                double newX = sumX / n_points;
                double newY = sumY / n_points;
                double newZ = sumZ / n_points;
                Point3D newCentroid=new Point3D(newX,newY,newZ);
                centroid.setPositionPoint(newCentroid);
            }
        }
    }
}