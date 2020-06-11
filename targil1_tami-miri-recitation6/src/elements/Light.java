package elements;

import primitives.Color;
/**
 *abstract class that all the classes that represent light source using it
 */
public abstract class Light {
	protected Color _intensity;
	/**
	 * constructor
	 * @param color _intensity value
	 */
	public Light(Color color)
	{
		_intensity=color;
	}
	/** 
	 * _intensity getter
	 * @return light color
	 */
	public Color getIntesity()
	{
		return _intensity;
	}

}
