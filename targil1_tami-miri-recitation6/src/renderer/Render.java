package renderer;


import java.util.List;

import elements.AmbientLight;
import elements.Camera;
import elements.LightSource;
import geometries.Intersectable;
import primitives.Color;
import primitives.Material;
import primitives.Point3D;
import primitives.Ray;
import primitives.Vector;
import scene.Scene;
import geometries.Intersectable.GeoPoint;
import static primitives.Util.*;
public class Render {
	
	private ImageWriter image;
	private Scene scene;
	/**
	 * ctr
	 * @param i image value
	 * @param s scene value
	 */
	public Render(ImageWriter i,Scene s)
	{
		image=i;
		scene=s;
	}
	/**
	 * Implements the image on the screen
	 */
	public void renderImage()
	{
		 Camera camera = scene.getCamera();
	     Intersectable geometries = scene.getGeometries();
	     java.awt.Color background = scene.getBackground().getColor();
	     AmbientLight ambientLight = scene.getAmbientLight();
	     double distance = scene.getDistance();
	     
	     int Nx = image.getNx();
	     int Ny = image.getNy();
	     double width = image.getWidth();
	     double height = image.getHeight();

	     for (int row = 0; row < Ny; row++) {
	         for (int collumn = 0; collumn < Nx; collumn++)
	         {
	                Ray ray = camera.constructRayThroughPixel(Nx, Ny, collumn, row, distance, width, height);
	                List<GeoPoint> intersectionPoints = geometries.findIntersections(ray);
	                if (intersectionPoints == null) {
	                    image.writePixel(collumn, row, background);
	                } else {
	                    GeoPoint closestPoint = getClosestPoint(intersectionPoints);
	                    image.writePixel(collumn, row, calcColor(closestPoint).getColor());
	                }
	         }
	     }
	}
	/**
	 * printing the grid
	 * @param interval interval between one square to an other
	 * @param color the grid color value
	 */
	public void printGrid(int interval, java.awt.Color color)
	{
		 int Nx = image.getNx();
	        int Ny = image.getNy();
	        for (int i = 0; i < Ny; i++) {
	            for (int j = 0; j < Nx; j++) {
	                if (i % interval == 0 || j % interval == 0) {
	                    image.writePixel(j, i, color);
	                }
	            }
	        }
	}
	/**
	 *call to WtiteToImage func in ImageWriter
	 */
	public void writeToImage() 
	{
	     image.writeToImage();
    }
	/**
	 * calculate the color in a given point
	 * @param p point value
	 * @return the intensity that defined in ambient lighting
	 */
	private Color calcColor(GeoPoint intersection)
	{
		Color result = scene.getAmbientLight().getIntesity();
        result = result.add(intersection.getGeometry().getEmmission());
        List<LightSource> lights = scene.getLights();

        Vector v = intersection.getPoint().subtract(scene.getCamera().getP0()).normalize();
        Vector n = intersection.getGeometry().getNormal(intersection.getPoint());

        Material material = intersection.getGeometry().getMaterial();
        int nShininess = material.getShininess();
        double kd = material.getKd();
        double ks = material.getKs();
        if (scene.getLights() != null) {
            for (LightSource lightSource : lights) {

                Vector l = lightSource.getL(intersection.getPoint());
                double nl = alignZero(n.dotProduct(l));
                double nv = alignZero(n.dotProduct(v));

                if (sign(nl) == sign(nv)) {
                    Color ip = lightSource.getIntensity(intersection.getPoint());
                    result = result.add(
                            calcDiffusive(kd, nl, ip),
                            calcSpecular(ks, l, n, nl, v, nShininess, ip)
                    );
                }
            }
        }

        return result;
    }
	
	private Object sign(double dotProduct) {
		return (dotProduct > 0d);
	}
	private Color calcDiffusive(double kd, double nl,Color ip) {
		if (nl < 0) nl = -nl;
        return ip.scale(nl * kd);
	}
	
	private Color calcSpecular(double ks,Vector l,Vector n, double nl, Vector v, int nShininess, Color lightIntensity) {
		
		double p = nShininess;
        Vector R = l.add(n.scale(-2 * nl)); // nl must not be zero!
        double minusVR = -alignZero(R.dotProduct(v));
        if (minusVR <= 0) {
            return Color.BLACK; // view from direction opposite to r vector
        }
        return lightIntensity.scale(ks * Math.pow(minusVR, p));
	}
	
	
	/**
	 * @param points list of all the points
	 * @return the closest cutting point with p0
	 */
	private GeoPoint getClosestPoint(List<GeoPoint> points)
	{
		Point3D p0= scene.getCamera().getP0();//p0 point
		double dis;
		double min=Double.MAX_VALUE;// the minimal distance
		GeoPoint p=null;
		for(GeoPoint geopoint : points)//for every point in the list we check if the distance from p0 is less then min
		{
			dis=p0.distance(geopoint.getPoint());
			if(dis<min)
			{
				min=dis;
				p=geopoint;
			}
		}
		return p;
	}
	

}
