package geometries;

import primitives.Ray;

import java.util.List;

import primitives.Point3D;


public interface Intersectable 
{
	/**
	 * @param ray
	 * @return a point3D list of itersections with a given ray
	 */
	public List<GeoPoint> findIntersections (Ray ray);
/**
 * class with Geometry and Point3D values with ctr, getter funcs and  equals func
 */
	public static class GeoPoint {
	    public Geometry geometry;
	    public Point3D point;
	   
	    /**
	      * ctr
	      * @param g geometry value
	      * @param p point value
	      */
	    public GeoPoint(Geometry g,Point3D p)
	    {
	    	geometry=g;
	    	point=p;
	    }
	    /**
	     * point getter
	     * @return point value
	     */
	    public Point3D getPoint() {
	    	return point;
	    }
	    /**
	     * geometry getter
	     * @return geomtry value
	     */
	    public Geometry getGeometry() {
	    	return geometry;
	    }
	    
	    @Override
	    /**
	     * func that check that the comparison is made on the same geometries and the points are the same
	     */
	    public boolean equals(Object obj) {
	       if (this == obj) return true;
	       if (obj == null) return false;
	       if (!(obj instanceof GeoPoint)) return false;
	       GeoPoint oth = (GeoPoint)obj;
	       return point.equals(oth.point)&& geometry.equals(oth.geometry);
	    }
	}


}
