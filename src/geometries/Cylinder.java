package geometries;

import java.util.List;

import primitives.Color;
import primitives.Material;
import primitives.Point3D;
import primitives.Ray;
import primitives.Vector;
/*
 * 
 */
public class Cylinder extends Tube
{
	double _height;
/**
 * ctr with a default _emmission value
 * @param r _radius value
 * @param ray axisRay value
 * @param height _height value
 */
	public Cylinder(double r, Ray ray,double height) {
		this(Color.BLACK,new Material(0, 0, 0),r,ray,height);
	}
	/**
	 * ctr
	 * @param emmission _emmission value
	 * @param r _radius value
     * @param ray axisRay value
     * @param height _height value
	 */
	public Cylinder(Color emmission,double r, Ray ray,double height) {
		this(emmission,new Material(0, 0, 0),r,ray,height);
	}
	/**
	 * calling Tube ctr
	 * @param emmission _emmission value
	 * @param material  _material
	 * @param r _radius value
	 * @param ray axisRay value
	 * @param height _hight value
	 */
	public Cylinder(Color emmission, Material material,double r, Ray ray,double height) {
		super(emmission,material,r, ray);
		// TODO Auto-generated constructor stub
		_height=height;
	}
	public Vector getNormal(Point3D _point) 
	{
	        return super.getNormal(_point);
	}
	/**
	 * _hight getter
	 * @return
	 */
	public double getH()
	{
		return _height;
	}
	@Override
	public String toString() 
	{
		 return "height:"+_height;
	} 
	    
	@Override
	public List<GeoPoint> findIntersections(Ray ray) {
		// TODO Auto-generated method stub
		return super.findIntersections(ray);
	}
	
	

}
