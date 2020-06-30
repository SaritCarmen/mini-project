package elements;

import primitives.Color;
import primitives.Point3D;
import primitives.Vector;
import static primitives.Util.*;

public class SpotLight extends PointLight {
	/**
	 * ctr that calls Light ctr with color value and vector
	 * @param color _intensity value
	 * @patam v direction value
	 * @param _position _position value
	 * @param kc _kC value 
	 * @param kl _kL value
	 * @param kq _kQ value
	 */
	private Vector _direction;
	public SpotLight(Color color,Vector v,Point3D _position,double _kC, double _kL, double _kQ,double r) {
		super(color, _position, _kC, _kL, _kQ,r);
		_direction=new Vector(v).normalize();
		// TODO Auto-generated constructor stub
	}
	/**
     * @return spotlight intensity
     */
    @Override
    public Color getIntensity(Point3D p) {
    	double projection = _direction.dotProduct(getL(p));

        if (isZero(projection)) {
            return Color.BLACK;
        }
        double factor = Math.max(0, projection);
        Color pointlightIntensity = super.getIntensity(p);
        return (pointlightIntensity.scale(factor));
    }
    

}
