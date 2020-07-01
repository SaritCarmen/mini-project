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
	 boolean IsIntersectionBox(Ray ray);
	 Box get_box() ;
	
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

	/**
	 * 
	 * class that united few geometries together in order to make the plan faster
	 *
	 */
	  public static class Box {

		  protected final double x0;
	        protected final double x1;
	        
	        protected final double y0;
	        protected final double y1;
	        
	        protected final double z0;
	        protected final double z1;

			/**
			 * @return the x0
			 */
			public double getX0() {
				return x0;
			}

			/**
			 * @return the x1
			 */
			public double getX1() {
				return x1;
			}

			/**
			 * @return the y0
			 */
			public double getY0() {
				return y0;
			}

			/**
			 * @return the y1
			 */
			public double getY1() {
				return y1;
			}

			/**
			 * @return the z0
			 */
			public double getZ0() {
				return z0;
			}

			/**
			 * @return the z1
			 */
			public double getZ1() {
				return z1;
			}

			
			/**
			 * @param x0
			 * @param x1
			 * @param y0
			 * @param y1
			 * @param z0
			 * @param z1
			 */
			public Box(double x0, double x1, double y0, double y1, double z0, double z1) {
				super();
				this.x0 = x0;
				this.x1 = x1;
				this.y0 = y0;
				this.y1 = y1;
				this.z0 = z0;
				this.z1 = z1;
			}
			
	        public boolean IntersectionBox(Ray r)
	        {
	        	double tmin = (x0 - r.getP().get_x().get()) / r.getV().get().get_x().get(); 
	        	double tmax = (x1 - r.getP().get_x().get()) / r.getV().get().get_x().get(); 
	    	 
	    	    if (tmin > tmax) {
	    	    	double temp=tmin;
	    			tmin=tmax;
	    			tmax=temp;
	    	    }
	    	    
	    	    double tymin = (y0 - r.getP().get_y().get()) / r.getV().get().get_y().get(); 
	        	double tymax = (y1 - r.getP().get_y().get()) / r.getV().get().get_y().get(); 
	    	 
	    	    if (tymin > tymax) {
	    	    	double temp=tymin;
	    			tymin=tymax;
	    			tymax=temp;
	    	    }    	
	    	 
	    	    if ((tmin > tymax) || (tymin > tmax)) 
	    	        return false; 
	    	 
	    	    if (tymin > tmin) 
	    	        tmin = tymin; 
	    	 
	    	    if (tymax < tmax) 
	    	        tmax = tymax; 
	    	 
	    	    double tzmin = (z0 - r.getP().get_z().get()) / r.getV().get().get_z().get(); 
	        	double tzmax = (z1 - r.getP().get_z().get()) / r.getV().get().get_z().get(); 
	    	 
	    	    if (tzmin > tzmax) {
	    	    	double temp=tzmin;
	    			tzmin=tzmax;
	    			tzmax=temp;
	    	    }
	    	 
	    	    if ((tmin > tzmax) || (tzmin > tmax)) 
	    	        return false; 
	    	 
	    	    if (tzmin > tmin) 
	    	        tmin = tzmin; 
	    	 
	    	    if (tzmax < tmax) 
	    	        tmax = tzmax; 
	    	 
	    	    return true;     	           
	        }

	        
	    }

}
