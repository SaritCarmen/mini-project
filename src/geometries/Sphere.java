package geometries;

import java.util.List;

import primitives.Color;
import primitives.Material;
import primitives.Point3D;
import primitives.Ray;
import primitives.Vector;
import static primitives.Util.*;

public class Sphere extends RadialGeometry
{
    Point3D _center;
	/**
	 * ctr with default emission and material values
	 * @param r _radius value
	 * @param _p _center value
	 */
	public Sphere(double r,Point3D _p)
	{
		this(Color.BLACK,new Material(0, 0, 0),r,_p);
    }
	/**
	 * ctr with default material value
	 * @param emmission _emission value
	 * @param r _radius value
	 * @param _p _center value
	 */
	public Sphere(Color emmission,double r,Point3D _p)
	{
		this(emmission,new Material(0,0,0),r,_p);
    }
	/**
	 * ctr that calls  RadialGeometry ctr with emission, material and radius values
	 * @param emmission _emisiion value
	 * @param material _material value
	 * @param r _radius value
	 * @param _p _center value
	 */
	public Sphere(Color emmission, Material material,double r,Point3D _p)//ctr
	{
		super(emmission,material,r);
		// TODO Auto-generated constructor stub
		_center=new Point3D(_p);
    }
/**
 * func that get a point and return the normal to the point
 */
	public Vector getNormal(Point3D _point) 
	{
	    Vector normal=_point.subtract(_center);
	    return normal.normalize();
	}
	//get point _center
	public Point3D getC()
	{
		return _center;
	}
	 @Override
	    public String toString() 
		{
		        return "center:"+_center;
		} 
	 @Override
	 /**
	  * @param ray
      * @return list of Point3D that intersect the given ray with the sphere
	  */
		public List<GeoPoint> findIntersections (Ray ray) {
			// TODO Auto-generated method stub
			  Point3D p0 = ray.getP();
		        Vector v = ray.getV();
		        Vector u;
		        try {
		            u = _center.subtract(p0);   // p0 == _center
		        } catch (IllegalArgumentException e) {
		            return List.of(new GeoPoint(this,ray.getTargetPoint(_radius)));
		        }
		        double tm = alignZero(v.dotProduct(u));
		        double dSquared = (tm == 0) ? u.lengthSquared() : u.lengthSquared() - tm * tm;
		        double thSquared = alignZero(_radius * _radius - dSquared);

		        if (thSquared <= 0) return null;

		        double th = alignZero(Math.sqrt(thSquared));
		        if (th == 0) return null;

		        double t1 = alignZero(tm - th);
		        double t2 = alignZero(tm + th);
		        if (t1 <= 0 && t2 <= 0) return null;
		        if (t1 > 0 && t2 > 0) return List.of(new GeoPoint(this,ray.getTargetPoint(t1)), new GeoPoint(this,ray.getTargetPoint(t2))); //P1 , P2
		        if (t1 > 0)
		            return List.of(new GeoPoint(this,ray.getTargetPoint(t1)));
		        else
		            return List.of(new GeoPoint(this, ray.getTargetPoint(t2)));
		    }
		}
	    

	


