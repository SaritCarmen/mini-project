package geometries;

import java.util.ArrayList;

import java.util.List;

import primitives.Ray;

public class Geometries implements Intersectable
{
     
	private List<Intersectable> _geometries = new ArrayList<Intersectable>();
	private Box _box;

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
        createBox();
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
    @Override
	public boolean IsIntersectionBox(Ray ray) {
		
		return this._box.IntersectionBox(ray);
	}
	@Override
	public Box get_box() {
		
		return _box;
	}
	  private void createBox() {
	        double x1=Double.NEGATIVE_INFINITY;
			double x0=Double.POSITIVE_INFINITY;
			double y1=Double.NEGATIVE_INFINITY;
			double y0=Double.POSITIVE_INFINITY;
			double z1=Double.NEGATIVE_INFINITY;
			double z0=Double.POSITIVE_INFINITY;
	        for(Intersectable geo: _geometries) { 
	        	
	           
	        	if(geo.get_box().getX0()<x0)
	        		x0=geo.get_box().getX0();
	        	if(geo.get_box().getX1()>x1) 
	        		x1=geo.get_box().getX1();
	        	if(geo.get_box().getY0()<y0) 
	        		y0=geo.get_box().getY0();
	        	if(geo.get_box().getY1()>y1) 
	        		y1=geo.get_box().getY1();
	        	if(geo.get_box().getZ0()<z0)
	        		z0=geo.get_box().getZ0();
	        	if(geo.get_box().getZ1()>z1) 
	        		z1=geo.get_box().getZ1();
	        	
	        }
	        
	        this._box=new Box(x0,x1,y0,y1,z0,z1);
	    }

}
