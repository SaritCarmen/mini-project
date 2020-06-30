package k_means;


import java.util.Random;

import geometries.*;
import primitives.*;

/**
 * @author rinat
 *
 */
public class Point {
	private Geometry geometry;
    private Point3D PositionPoint;
    private int cluster_number = 0;

    public Point(Geometry g)
    {
        this.geometry=g;
        this.PositionPoint=g.getPositionPoint();
    }
   

    /**
	 * @return the geometry
	 */
	public Geometry getGeometry() {
		return geometry;
	}


	/**
	 * @param geometry the geometry to set
	 */
	public void setGeometry(Geometry geometry) {
		this.geometry = geometry;
	}


	/**
	 * @return the positionPoint
	 */
	public Point3D getPositionPoint() {
		return PositionPoint;
	}


	/**
	 * @param positionPoint the positionPoint to set
	 */
	public void setPositionPoint(Point3D positionPoint) {
		PositionPoint = positionPoint;
	}


	/**
	 * @return the cluster_number
	 */
	public int getCluster_number() {
		return cluster_number;
	}


	/**
	 * @param cluster_number the cluster_number to set
	 */
	public void setCluster_number(int cluster_number) {
		this.cluster_number = cluster_number;
	}


	//Calculates the distance between two points.
    protected static double distance(Point p, Point centroid) {
    	return p.getPositionPoint().distanceSquared(centroid.getPositionPoint());
    }

    //Creates random point
    protected static Point createRandomPoint() {
    	return Point.createRandomPoint(Double.NEGATIVE_INFINITY,Double.POSITIVE_INFINITY);
    }
    
    protected static Point createRandomPoint(double min, double max) {
        Random r = new Random();
        double x = min + (max - min) * r.nextDouble();
        double y = min + (max - min) * r.nextDouble();
        double z = min + (max - min) * r.nextDouble();
        return new Point(new Sphere(1,new Point3D(x,y,z)));
    }
    /*
    protected static List createRandomPoints(int min, int max, int number) {
        List points = new ArrayList(number);
        for(int i = 0; i < number; i++) {
            points.add(createRandomPoint(min,max));
        }
        return points;
    }*/

    public String toString() {
        return "("+this.PositionPoint.get_x()+","+this.PositionPoint.get_y()+","+this.PositionPoint.get_z()+")";
    }
}