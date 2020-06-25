package primitives;
/**
 *class that have to every geomtry which represent the rank  of the brightness of the geometry
 *
 */
public class Material {
	private double _kD;//distance factor
	private double _kS;//distance factor
	private int _nShininess;//the rank  of the brightness of the material
	private double _kT;//Promotes transparency
	private double _kR;//Reflection coefficient.
	
	
	/**
	 *Material constructor
	 * @param kd _kD value
	 * @param ks _kS value
	 * @param shininess _nShininess value
	 */
	public Material(double kd, double ks, int shininess) {
		this(kd,ks,shininess,0,0);
	}
	
	public Material(double kd, double ks, int shininess,double kT,double kR) {
		_kD=kd;
		_kS=ks;
		_nShininess=shininess;
		_kR =kR;
		_kT =kT;
	}
	
	/** 
	 * @return_kD value
	 */
	public double getKd() {
		return +_kD;
	}
	/**
	 * @return _kS value
	 */
	public double getKs() {
		return _kS;
	}
	/**
	 * @return _nShininess value
	 */
	public int getShininess() {
		return _nShininess;
	}
	
	public double getkR() {
		return _kR;
	}
	
	public double getkT() {
		return _kT;
	}
}
