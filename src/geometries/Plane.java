package geometries;


import primitives.Color;
import primitives.Material;

import java.util.List;
import static primitives.Util.*;

import primitives.Point3D;
import primitives.Ray;
import primitives.Vector;

public class Plane extends Geometry
{
	Point3D p1;
	Vector normal;
	/**
	 * Plane  constructor receiving 3 Point3D value
     * @param p1  point p1
	 * @param p2 vector normal
	 * @param p3 vector normal
	 */
	public Plane(Point3D p1,Point3D p2, Point3D p3)
	{
		super();
		this.p1=p1;//relation point of the plane
		Vector v1=p2.subtract(p1);
		Vector v2=p3.subtract(p2);
		normal=(v1.crossProduct(v2)).normalize();
		
        
	}
	/**
	 * Plane  constructor receiving Point3D value and Vector
	 * @param _p point p1
	 * @param v Vector normal
	 */	
	public Plane(Point3D _p,Vector v)
	{
		super();
		p1=new Point3D(_p);
		normal=new Vector(v);
	}
	/**
	 * ctr that calls Geometry ctr with emission and material valus 
	 * @param emmission _emission value
	 * @param material _material value
	 * @param p1 point p1
	 * @param p2 vector normal
	 * @param p3 vector normal
	 */
	public Plane(Color emmission, Material material,Point3D p1,Point3D p2, Point3D p3 )
	{
		super(emmission, material);
		this.p1=p1;//relation point of the plane
		Vector v1=p2.subtract(p1);
		Vector v2=p3.subtract(p2);
		normal=(v1.crossProduct(v2)).normalize();

	}
	/**
	 * ctr with default material
	 * @param emmission _emission value
	 * @param p1 point p1
	 * @param p2  vector normal
	 * @param p3 vector normal
	 */
	public Plane(Color emmission,Point3D p1,Point3D p2, Point3D p3 )
	{
		this(emmission ,new Material(0,0,0),p1,p2,p3);
	}
	
    public Point3D getP()//get point p1
    {
    	return p1;
    }
    
    public Vector getV()//get vector normal
    {
    	return normal;
    }
    
    @Override
    public String toString() 
	{
	        return "point:"+p1.toString()+"normal vector:"+normal.toString();
	}

	@Override
	public Vector getNormal(Point3D p) {
		// TODO Auto-generated method stub
		return this.normal;
	}
	@Override
	/**
	 * @param ray
     * @return list of Point3D that intersect the given ray with the plane
	 */
	public List<GeoPoint> findIntersections (Ray ray)
	{
		// TODO Auto-generated method stub
		  Vector p0Q;
	        try {
	            p0Q = p1.subtract(ray.getP());
	        } catch (IllegalArgumentException e) {
	            return null; // ray starts from point Q - no intersections
	        }

	        double nv = normal.dotProduct(ray.getV());
	        if (isZero(nv)) // ray is parallel to the plane - no intersections
	            return null;

	        double t = alignZero(normal.dotProduct(p0Q) / nv);

	        return t <= 0 ? null : List.of(new GeoPoint(this,ray.getTargetPoint(t)));
	 }



}
