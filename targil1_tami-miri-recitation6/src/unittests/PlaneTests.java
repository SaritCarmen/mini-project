/**
 * 
 */
package unittests;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

import geometries.Intersectable.GeoPoint;
import geometries.Plane;
import primitives.Point3D;
import primitives.Ray;
import primitives.Vector;

/**
 * testing getNormal() of plane
 *
 */
public class PlaneTests {

	/**
	 * Test method for {@link geometries.Plane#getNormal(primitives.Point3D)}.
	 */
	@Test
	public void testGetNormal() {
		 // ============ Equivalence Partitions Tests ==============
		 Plane pl = new Plane(new Point3D(0, 0, 1), new Point3D(1, 0, 0), new Point3D(0, 1, 0));
	     double sqrt3 = Math.sqrt(1d / 3);
	     assertEquals("Bad normal to plane", new Vector(sqrt3, sqrt3, sqrt3), pl.getNormal(new Point3D(0, 0, 1)));
	}
	
	@Test
	public void testFindIntersections() 
	{
		Plane plane =new Plane(new Point3D(1,1,1), new Vector(1,0,0));

		// ============ Equivalence Partitions Tests ==============

		        // TC01: Ray intersects the plane (1 points)

		        GeoPoint p1 = new GeoPoint(plane,new Point3D(1,2,2));
		        List<GeoPoint> result = plane.findIntersections(new Ray(new Point3D(2,2,2),new Vector(-0.5,0,0)));
		        assertEquals("Wrong number of points",1, result.size());
		        assertEquals("Ray crosses plane", List.of(p1), result);
		       
		     // TC02: Ray does not intersect the plane
		       
		        List<GeoPoint> result1 = plane.findIntersections(new Ray(new Point3D(2,2,2),new Vector(0.5,0,0)));
                assertEquals("Ray crosses plane", null, result1);
		       
		     // =============== Boundary Values Tests ==================
		       
		     // TC03:Ray is parallel to the plane,the ray not included in the plan
		       
		        List<GeoPoint> result2 = plane.findIntersections(new Ray(new Point3D(2,2,2),new Vector(0,1,1)));
		        assertEquals("Ray crosses plane", null, result2);
		    
		     // TC04: Ray is orthogonal to the plane and before the plane

		        List<GeoPoint> result3 = plane.findIntersections(new Ray(new Point3D(2,2,2),new Vector(1,0,0)));
		        assertEquals("Ray crosses plane", null, result3);
		   
		      // TC05: Ray is orthogonal to the plane and in the plane

		        List<GeoPoint> result4 = plane.findIntersections(new Ray(new Point3D(2,2,2),new Vector(-1,0,0)));
		        assertEquals("Wrong number of points",1, result4.size());
		        assertEquals("Ray crosses plane", List.of(new GeoPoint(plane,new Point3D(1.0,2.0,2.0))), result4);
		       
		      // TC06: Ray is orthogonal to the plane and after the plane

		        List<GeoPoint> result5 = plane.findIntersections(new Ray(new Point3D(-2,-2,-2),new Vector(-1,0,0)));
		        assertEquals("Ray crosses plane", null, result5);
		       
		     // TC07: Ray is neither orthogonal nor parallel to and begins at the plane
		       
		        List<GeoPoint> result6 = plane.findIntersections(new Ray(new Point3D(1,2,2),new Vector(1,0,1)));
		        assertEquals("Ray crosses plane", null, result6);
		       
		       
		
	}

}
