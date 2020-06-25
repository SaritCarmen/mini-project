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
	 * Fixed for moving first rays for shading, transparency and reflection rays 
	 * its value can be run instead of orders of magnitude.
	 */
	private static final double DELTA = 0.1;
	/**
	 * stop condition  maximum level for the recursion
	 */
	private static final int MAX_CALC_COLOR_LEVEL = 10;
	/**
	 *  stop condition  minimum level for the recursion
	 */
	private static final double MIN_CALC_COLOR_K = 0.001;
	
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
	                GeoPoint closestPoint = findClosestIntersection(ray);

	               // List<GeoPoint> closestPoint = geometries.findIntersections(ray);
	                if (closestPoint == null) {
	                    image.writePixel(collumn, row, background);
	                } else {
	                   //  closestPoint = getClosestPoint(intersectionPoints);
	                    image.writePixel(collumn, row, calcColor(closestPoint,ray).getColor());
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
	 * @param geopoint point value
	 * @param inRay Ray value
	 * @return send to another function - recorsya
	 */
	private Color calcColor(GeoPoint geopoint, Ray inRay) {
		return calcColor(geopoint, inRay, MAX_CALC_COLOR_LEVEL, 1.0).add(
		scene.getAmbientLight().getIntesity());
		}

	
	/**
	 * calculate the color in a given point
	 * @param p point value
	 * @return the intensity that defined in ambient lighting
	 */
	private Color calcColor(GeoPoint intersection, Ray inRay, int level, double k )
	{
		if (level == 0 || k < MIN_CALC_COLOR_K)
			return Color.BLACK;
		Color color = intersection.geometry.getEmmission(); // remove Ambient Light
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
                	double ktr = transparency(lightSource, l, n, intersection);
                	if (ktr * k > MIN_CALC_COLOR_K) {
                		Color ip = lightSource.getIntensity(intersection.getPoint()).scale(ktr);;
                		color = color.add(
                            calcDiffusive(kd, nl, ip),
                            calcSpecular(ks, l, n, nl, v, nShininess, ip)                	
                				);
                	}
                }
            }
        }

        if (level == 1)
        	return Color.BLACK;
        double kr = intersection.geometry.getMaterial().getkR(), kkr = k * kr;
        if (kkr > MIN_CALC_COLOR_K) {
        	Ray reflectedRay = constructReflectedRay(n, intersection.point, inRay);
        	GeoPoint reflectedPoint = findClosestIntersection(reflectedRay);
        	if (reflectedPoint != null)
        		color = color.add(calcColor(reflectedPoint, reflectedRay,
        		level-1, kkr).scale(kr));
        }
        double kt = intersection.geometry.getMaterial().getkT();
        double kkt = k * kt;
        if (kkt > MIN_CALC_COLOR_K) {
        	Ray refractedRay = constructRefractedRay(intersection.point, inRay,n) ;
        	GeoPoint refractedPoint = findClosestIntersection(refractedRay);
        	if (refractedPoint != null)
        		color = color.add(calcColor(refractedPoint, refractedRay,level-1, kkt).scale(kt));
        }
        
        return color;
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
	
	/**
	 * Non-shading test between point and light source
	 * @param light the current light source
	 * @param l vector from light source
	 * @param n normal to raise the point in £ to fix the floating point problem
	 * @param gp point on the geometry which the vector cut
	 * @return true if the point is not hiding and false if it is
	 */
	private boolean unshaded(LightSource light,Vector l, Vector n, GeoPoint gp)
	{
		if(gp.getGeometry().getMaterial().getkT() != 0)
			return false;
		Vector lightDirection = l.scale(-1); // from point to light source
		Vector delta = n.scale(n.dotProduct(lightDirection) > 0 ? DELTA : - DELTA);
		Point3D point = gp.point.add(delta);
		Ray lightRay = new Ray(point, lightDirection);
		List<GeoPoint> intersections = scene.getGeometries().findIntersections(lightRay);
		if(intersections ==null)
			return true;
		double lightDistance = light.getDistance(point);
		for (GeoPoint geo : intersections) {
			 if (alignZero( geo.point.distance(point)-lightDistance) <= 0 )
				   return false;
		}
		return true;

	}

	/**
	 * Partial shading in case the body or bodies that block the light source from the point have transparency at some level	 * @param ls
	 * @param l
	 * @param n
	 * @param geopoint
	 * @return
	 */
	private double transparency(LightSource light, Vector l, Vector n, GeoPoint geopoint) {
		Vector lightDirection = l.scale(-1); // from point to light source
		Ray lightRay = new Ray(geopoint.getPoint(), lightDirection);
		List<GeoPoint> intersections = scene.getGeometries().findIntersections(lightRay);
		if(geopoint.getGeometry().getMaterial().getkT() != 0)
			return 1.0;
		Vector delta = n.scale(n.dotProduct(lightDirection) > 0 ? DELTA : - DELTA);
		Point3D point = geopoint.point.add(delta);
		if(intersections ==null)
			return 1.0;
		double ktr = 1.0;
		double lightDistance = light.getDistance(geopoint.getPoint());
		for (GeoPoint geo : intersections) {
			 if (alignZero( geo.point.distance(point)-lightDistance) <= 0 )
			 {
				 ktr *= geopoint.geometry.getMaterial().getkT();
				 if (ktr < MIN_CALC_COLOR_K) 
					 return 0.0;
			 }
		}
		return ktr;
		
	}
	
	/**
	 * find the closest intersection of the Ray and geometry.
	 * @param ray
	 * @return closest intersection of the Ray and geometry.
	 */
	private GeoPoint findClosestIntersection(Ray ray) {
		if (ray == null) 
		{
            return null;
        }
		Point3D point = ray.getP();
		List<GeoPoint> intersections = scene.getGeometries().findIntersections(ray);
		if( intersections == null)
			return null;
		double dis;
		double closestDistance =  Double.MAX_VALUE;
		GeoPoint p = null;
		for(GeoPoint geopoint : intersections)//for every point in the list we check if the distance from p0 is less then min
		{
			dis=point.distance(geopoint.getPoint());
			if(dis<closestDistance)
			{
				closestDistance=dis;
				p=geopoint;
			}
		}
		return p;
		
	}

	/**
	 * calculate reflection ray (R=D-2(D*N)N)
	 * @param pointGeo p value in ray ctr
	 * @param inRay source light direction vector
	 * @param n the mormal to the point
	 * @return ray after moving the start point of the ray on the normal to the geometry in the direction of the new ray	
*/
	private Ray constructReflectedRay(Vector n, Point3D point, Ray inRay) {
		Vector v = inRay.getV();
		double dotProduct = (v.dotProduct(n));
		if (dotProduct == 0) {
            return null;
        }
		Vector ans = v.subtract(n.scale(2*dotProduct));
		return new Ray(point,ans,n); 
	}

	private Ray constructRefractedRay( Point3D point, Ray inRay,Vector n) {
		Vector dir = inRay.getV();
		return new Ray(point,dir,n);
	}
	
	
 }
