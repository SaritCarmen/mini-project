package geometries;

import primitives.Color;
import primitives.Material;

public abstract class RadialGeometry extends Geometry
{
	double _radius;
	/** 
     * @return _radius value
     */
	public double getR() {
		return _radius;
	}
	/**
	 *ctr 
	 * @param r _radius value
	 */
	public RadialGeometry (double r)
	{
		_radius=r;
	}
	/**
	 * ctr that calls Geometry ctr with emission value
	 * @param emmission _emission value
	 * @param r _radius value
	 */
	public RadialGeometry (Color emmission,double r)
	{
		super(emmission);
		_radius=r;
	}
	/**
	 * ctr that calls Geometry ctr with emission and material values
	 * @param emmission _emission value
	 * @param material _material value
	 * @param r _radius value
	 */
	public RadialGeometry (Color emmission,Material material,double r)
	{
		super(emmission,material);
		_radius=r;
	}
	/**
	 * coppy ctr
	 * @param r _radius value
	 */
	public RadialGeometry (RadialGeometry  r)
	{
		_radius=r.getR();
	}

	

}
