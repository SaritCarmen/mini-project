package unittests;

import static org.junit.Assert.*;


import org.junit.Test;

import geometries.Tube;
import primitives.Point3D;
import primitives.Ray;
import primitives.Vector;
/**
 * testing getNormal() of tube
 */
public class TubeTests {

	@Test
	public void testGetNormal() {
		 // ============ Equivalence Partitions Tests ==============
		Tube t = new Tube(1,new Ray(new Point3D(1.0,0.0,0.0), new Vector(0,1,0)));
        assertEquals("Wrong normal to Cylinder", new Vector(1.0, 0.0, 0.0), t.getNormal(new Point3D(2, 0, 0)));
	}
	
	@Test
	public void testFindIntersections() 
	{
		
	}

}
