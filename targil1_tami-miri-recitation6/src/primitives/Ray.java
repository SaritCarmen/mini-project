package primitives;
import static primitives.Util.*;

public class Ray
{
	/**
	 * The set of points on the straight that are on one side relative to a given point on the straight.
	 *  Defined by Vector and Point3D
	 */
	
	private Point3D p;
	private Vector dir;
	/**
	    * Ray constructor receiving point3D value and Vector
	    *  check it is normalized"
		 * 
		 * @param _p Point3D p
		 * @param v  Vector dir
		 * 
	   */
	
	public Ray(Point3D _p, Vector v)
	{
		p=new Point3D(_p);
		dir=new Vector(v).normalized();	
	}
	/**
	    * Ray  constructor receiving Ray value 
		 * 
		 * @param r Ray for constructor 
		 * 
	   */
	public Ray(Ray r)
	{
		p=new Point3D(r.p);
		dir=new Vector(r.dir).normalized();	
	}
	 /**
     * Point3D value getter
     * 
     * @return Point3D p 
     */
	public Point3D getP()
	{
		return p;
	}
	 /**
     * Vector value getter
     * 
     * @return Vector dir
     */
	public Vector getV()
	{
		return dir;
	}
	public boolean equals(Object obj) {
	      if (this == obj) return true;
	      if (obj == null) return false;
	      if (!(obj instanceof Ray)) return false;
	      Ray oth = (Ray)obj;
	      return p.equals(oth.p) && dir.equals(oth.dir);
	   }
	/**
	 * calculate point on a ray (p=p0+t*v)
	 * @param length as t value
	 * @return point3D on the ray
	 */
    public Point3D getTargetPoint(double length) 
    {
         return isZero(length ) ?p : p.add(dir.scale(length));
    }
	@Override
	public String toString() 
	{
	        return "point:"+p.toString()+"diraction vector:"+dir.toString();
	}
	
	
	


}
