package scene;

/**
*@author sarit carmen 315154294
*@author yael vaknin 207476268
*/

import elements.AmbientLight;
import elements.Camera;
import geometries.Geometries;
import geometries.Intersectable;
import primitives.Color;

import java.util.LinkedList;
import java.util.List;
import elements.*;

public class Scene 
{
	private String _name;
	private Color _background;
	private AmbientLight _ambientLight;
	private Geometries _geometries;
	private Camera _camera;
	private double _distance;
	private List<LightSource>_lights;
	
	/**
	 * ctr that create an empty collection of elements in 3D model. 
	 * @param name scene name
	 */
	public Scene(String name)
	{
		_name=name;
		_geometries=new Geometries();
	}
	/** 
	 * @return list of lights source
	 */
	public List<LightSource> getLights(){
		return _lights;
	}
	/** 
	 * @return scene name
	 */
	public String getName() {
		 return _name;
	}
	/**
	 * @return background color
	 */
	public Color getBackground()
	{
		return _background;
	}
	/**
	 * @return ambient light
	 */
	public AmbientLight getAmbientLight()
	{
		return _ambientLight;
	}
	/**
	 * @return camera
	 */
	public Camera getCamera()
	{
		return _camera;
	}
	/**
	 * @return distance from the camera to the screen
	 */
	public double getDistance()
	{
		return _distance;
	}
	/**
	 * @return collection of elements in 3D model.
	 */
	public Geometries getGeometries()
	{
		return _geometries;
	}
	/**
	 * @param c camera value
	 */
	public void setCamera(Camera c)
	{
		_camera=c;
	}
	/**
	 * @param d distance value
	 */
	public void setDistance(double d)
	{
		_distance=d;
	}
	/**
	 * @param b background color value
	 */
	public void setBackground(Color b)
	{
		_background=b;
	}
	/**
	 * @param a ambient light value
	 */
	public void setAmbientLight(AmbientLight a)
	{
		_ambientLight=a;
	}
	/**
	 * func that add geometries to the 3D model
	 * @param geometries geometries to add the collection 
	 */
	public void addGeometries(Intersectable...  geometries) 
	{
		for(Intersectable geo :geometries)
		{
			_geometries.add(geo);
		}
	}
	/**
	 * func that add light sorces to the list 
	 * @param lights list of light sources to add 
 	 */
	public void addLights(LightSource... lights) 
	{
		if(_lights==null)
			_lights=new LinkedList<LightSource>();
		
		for(LightSource light :lights)
		{
			_lights.add(light);
		}
		
	}
	

    


}
