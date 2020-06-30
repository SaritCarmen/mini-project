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
   protected Box _box;
   abstract public Point3D getPositionPoint();
   
   
   /**
    * Geomtry ctr with default _material value
    * @param color _emmission value
    */
   public Geometry(Color color,Box box)
   {
	   this (color,new Material(0, 0, 0) ,box);
	   
   }
   /**
    * default ctr
    */
   public Geometry(Box box) {
	   this(Color.BLACK,new Material(0, 0, 0),box );
	  }
   /**
    * ctr
    * @param color _emission value
    * @param material _material value
    */
   public Geometry(Color color, Material material,Box box)
   {
	   _emmission=color;
	   _material=material;
	   _box =box;
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
   
   /**
    * getter
    * @return _box
    */
   public Box get_box()
   {
	   return _box;
   }

/**
 * setter
 * @param box is _box
 */
   public void set_box(Box box)
   {
	   _box=box;
   }
   
}
