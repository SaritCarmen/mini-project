package elements;

import primitives.Color;

public class AmbientLight extends Light {
	
	/**
	 * ctr that calls Light ctr
	 * @param ia filling Light Intensity
	 * @param ka coefficient of attenuation of filler light
	 */
	public AmbientLight(Color ia, double ka)
	{
		super(ia.scale(ka));

	}

}
