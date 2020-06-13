package geometries;

import primitives.Color;
import primitives.Material;
import primitives.Point3D;
import static primitives.Util.*;
import java.util.List;

import primitives.Ray;
import primitives.Vector;

public class Tube extends RadialGeometry
{
    Ray axisRay;
	/**
	 * ctr with default _emmission value
	 * @param r _radius value
	 * @param ray axisRay value
	 */
	public Tube(double r,Ray ray) 
	{
		this(Color.BLACK,new Material(0, 0, 0),r,ray);
	}
	/**
	 * ctr
	 * @param emmission _emmission value
	 * @param r _radius value
	 * @param ray axisRay value
	 */
	public Tube(Color emmission,double r,Ray ray) 
	{
		this(emmission,new Material(0, 0, 0),r,ray);
	}
	/**
	 * calling RadialGeometry ctr
	 * @param emmission _emmission value
	 * @param material  _material
	 * @param r _radius value
	 * @param ray axisRay value
	 */
	public Tube(Color emmission,Material material,double r,Ray ray) 
	{
		super(emmission,material,r);
		// TODO Auto-generated constructor stub
		axisRay=new Ray(ray);
	}
	//get ray axisRay
	public Ray getRay()
	{
		return axisRay;
	}

	public Vector getNormal(Point3D _point) 
	{
	   Point3D o=axisRay.getP();
	   Vector v=axisRay.getV();
	   double t=_point.subtract(o).dotProduct(v);
	   if(!isZero(t))
	   {
		   o=o.add(v.scale(t));
	   }
	   return _point.subtract(o).normalize();
	   
	   
	 
	}
	
	
	public List<GeoPoint> findIntersections (Ray ray) {
		// TODO Auto-generated method stub
	
		    return null;
		
	} 
	
	@Override
	public String toString() 
    {
		 return "ray:"+axisRay;
	} 
	

}
