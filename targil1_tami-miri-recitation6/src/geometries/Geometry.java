package geometries;

import primitives.Color;
import primitives.Material;
import primitives.Point3D;
import primitives.Vector;

public abstract class Geometry implements Intersectable 
{
	/**
	 * func that calculate the normal vector to the object in the given point p
	 * @param p point3D value
	 * @return normal vector 
	 */
   abstract	public Vector getNormal(Point3D p);
   protected Color _emmission;
   protected Material _material ;
   /**
    * Geomtry ctr with default _material value
    * @param color _emmission value
    */
   public Geometry(Color color)
   {
	   _material=new Material(0, 0, 0) ;
	   _emmission=color;
   }
   /**
    * default ctr
    */
   public Geometry() {
	   _material=new Material(0, 0, 0) ;
	   _emmission=Color.BLACK;
   }
   /**
    * ctr
    * @param color _emission value
    * @param material _material value
    */
   public Geometry(Color color, Material material)
   {
	   _emmission=color;
	   _material=material;
   }
  /** 
   * emisiion getter
   * @return  _emmission
   */
   public Color getEmmission()
   {
	   return _emmission;
   }
   /**
    * _material getter
    * @return _material
    */
   public Material getMaterial() {
	   return _material;
   }
   
}
