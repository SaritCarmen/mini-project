package primitives;

public class Vector 
{
	/**
	 * Vector in the plane Defined by the end point (when starting point - First deputies).
	 */
	Point3D p;
   /**
    * Vector  constructor receiving 3 coordinates value and check it is not "the 0 vector"
	 * 
	 * @param a coordinate x value in point p
	 * @param b coordinate y value in point p
	 * @param c coordinate z value in point p
	 * 
   */
	public Vector(Coordinate a,Coordinate b,Coordinate c)throws IllegalArgumentException 
	{
		if(Point3D.ZERO.equals(new Point3D(a,b,c)))
			throw new IllegalArgumentException("the zero vector");
		p=new Point3D(a,b,c);
		
	}
	/**
	    * Vector  constructor receiving 3 double value and check it is not "the 0 vector"
		 * 
		 * @param a coordinate x value in point p
		 * @param b coordinate y value in point p
		 * @param c coordinate z value in point p
		 * 
	   */
	public Vector(double a,double b,double c)throws IllegalArgumentException 
	{
		if(Point3D.ZERO.equals(new Point3D(a,b,c)))
			throw new IllegalArgumentException("the zero vector");
		p=new Point3D(a,b,c);
	}
	/**
	    * Vector constructor receiving a point3D value and check it is not "the 0 vector"
		 * 
		 * @param _p p value
		 * 
	   */
	public Vector(Point3D _p)throws IllegalArgumentException 
	{
		if(_p.equals(Point3D.ZERO))
			throw new IllegalArgumentException("the zero vector");
		p= new Point3D(_p);
	}
	/**
	    * Vector  constructor receiving a vector value 
		 * @param v for v.p point value
	   */
	public Vector(Vector v)
	{
		p=new Point3D(v.p);
	}
	 /**
     * Point3D value getter
     * 
     * @return Point3D  value
     */
    public Point3D get() {
        return p;
    }
    
    @Override
    public boolean equals(Object obj) {
       if (this == obj) return true;
       if (obj == null) return false;
       if (!(obj instanceof Vector)) return false;
       Vector oth = (Vector)obj;
       return p.equals(oth.p);
    }
    @Override
    public String toString() {
        return "point:"+p.toString();
    }
    /**
     * add/subtract vectors
     * @param v second vector
     * @return new vector
     */

	public Vector subtract(Vector v)
	{
		Vector nv=new Vector(this.p.subtract(v.p));
		return nv;
	}
	public Vector add(Vector v)
	{
		Vector nv=new Vector(this.p.add(v));
		return nv;
	}
	/**
	 * double vector in scalar
	 * @param tmp scalar
	 * @return new vector
	 */	
	public Vector scale(double tmp)
	{
		return new Vector(p.x.get()*tmp,p.y.get()*tmp,p.z.get()*tmp);
	}
	/**
	 * dot-product
	 * @param vector dots values
	 * @return scalar
	 */
	 public double dotProduct(Vector vector) 
	 {
	        return p.get_x().get() * vector.p.get_x().get() +p.get_y().get() * vector.p.get_y().get() + p.get_z().get() * vector.p.get_z().get();
	 }
	 /**
	  * cross-product
	  * @param v dots value
	  * @return Vector that stands for the two existing vectors
	  */
	 public Vector crossProduct(Vector v)
	 {
		return new Vector(this.p.get_y().get()*v.p.get_z().get()-this.p.get_z().get()*v.p.get_y().get(),
				this.p.get_z().get()*v.p.get_x().get()-this.p.get_x().get()*v.p.get_z().get(),
				this.p.get_x().get()*v.p.get_y().get()-this.p.get_y().get()*v.p.get_x().get());
		 
	 }
	 /**
	  *@return Square vector length
	  */
	 public double lengthSquared()
	 {
		 return p.get_x().get()*p.get_x().get()+p.get_y().get()*p.get_y().get()+p.get_z().get()*p.get_z().get();
	 }
	 /**
	  * @return vector length
	  */
	 public double length()
	 {
		return Math.sqrt(this.lengthSquared()); 
	 }
	 /**
	  * Normalize the vector that will change the vector itself
	  * @return the vector
	  */
	 public Vector normalize()
	 {
		 this.p = new Point3D(this.scale(1/this.length()).get());
		 return this;
	 }
	 /**
	  * Normalize a vector 
	  * @return new vector normalized in the same direction as a vector
	  */
	 public Vector normalized()
	 {
		 Vector v=new Vector(this);
		 return v.normalize();
	 }
}
