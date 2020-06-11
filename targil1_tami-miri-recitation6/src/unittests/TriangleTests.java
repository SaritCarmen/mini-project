/**
 * 
 */
package unittests;

import static org.junit.Assert.*;
import geometries.Intersectable.GeoPoint;

import java.util.List;

import org.junit.Test;

import geometries.Polygon;
import geometries.Triangle;
import primitives.Point3D;
import primitives.Ray;
import primitives.Vector;

/**
 * testing getNormal() of triangle
 *
 */
public class TriangleTests {

	/**
	 * Test method for {@link geometries.Triangle#getNormal(primitives.Point3D)}.
	 */
	@Test
	public void testGetNormal() {
		 // ============ Equivalence Partitions Tests ==============
		Polygon p = new Polygon(new Point3D(0, 0, 1), new Point3D(1, 0, 0), new Point3D(0, 1, 0));
	     double sqrt3 = Math.sqrt(1d / 3);
	     assertEquals("Bad normal to triangle", new Vector(sqrt3, sqrt3, sqrt3), p.getNormal(new Point3D(0, 0, 1)));
	}
	
	@Test
	public void testFindIntersections() 
	{
		Triangle triangle =new Triangle(new Point3D(0,0,1),new Point3D(1,0,0), new Point3D(-1,0,0));
		 // ============ Equivalence Partitions Tests ==============
		
		 //TC01 Inside triangle
		 assertEquals("Ray Inside the triangle",List.of(new GeoPoint(triangle, new Point3D(0,0,0.5))),triangle.findIntersections((new Ray(new Point3D(0, 2, 0.5), new Vector(0, -1, 0)))));
		 //TC02 Outside against edge
		 assertEquals("Ray starts outside against edge",null,triangle.findIntersections((new Ray(new Point3D(0.5, -2, -1), new Vector(0, 1, 0)))));
		 //TC03 Outside against vertex
		 assertEquals("Ray starts outside against vertex",null,triangle.findIntersections((new Ray(new Point3D(1.5, -2, -0.2), new Vector(0, 1, 0)))));
		// =============== Boundary Values Tests ==================
		//TC04 the point is on edge
		 assertEquals("Ray's point is on the edge",null,triangle.findIntersections((new Ray(new Point3D(0.5, -2, 0), new Vector(0, 1, 0)))));
		//TC05 the point is in vertex
		 assertEquals("Ray's point is in vertex",null,triangle.findIntersections((new Ray(new Point3D(1, -1, 0), new Vector(0, 1, 0)))));
		//TC06 the point is On edge's continuation
		 assertEquals("Ray's point is On edge's continuation",null,triangle.findIntersections((new Ray(new Point3D(2, -2, 0), new Vector(0, 1, 0)))));
    }
	

}
