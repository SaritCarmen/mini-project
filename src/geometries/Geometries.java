package geometries;

import java.util.ArrayList;

import java.util.List;

import primitives.Ray;

public class Geometries implements Intersectable
{
     
	private List<Intersectable> _geometries = new ArrayList<Intersectable>();
	

    public Geometries(Intersectable... _geometries)//ctr
    {
        add( _geometries);
    }
    /**
     * add geometries to the collection
     * @param geometries
     */
        
    public void add(Intersectable... geometries) {
        for (Intersectable geo : geometries ) {
            _geometries.add(geo);
        }
    }

    /**
     * @param ray
     * @return list of Point3D that intersect the collection
     */
    @Override
    public List<GeoPoint> findIntersections(Ray ray) {
    	List<GeoPoint> intersections = null;

        for (Intersectable geo : _geometries) {
        	List<GeoPoint> tempIntersections = geo.findIntersections(ray);
            if (tempIntersections != null) {
                if (intersections == null)
                    intersections = new ArrayList<GeoPoint>();
                intersections.addAll(tempIntersections);
            }
        }
        return intersections;

     }
	

}
