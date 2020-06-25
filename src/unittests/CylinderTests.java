/**
 * 
 */
package unittests;

import static org.junit.Assert.*;



import org.junit.Test;

import geometries.Cylinder;
import primitives.Point3D;
import primitives.Ray;
import primitives.Vector;

/**
*@author sarit carmen 315154294
*@author yael vaknin 207476268
 * testing getNormal() of Cylinder
 */
public class CylinderTests {

	/**
	 * Test method for {@link geometries.Cylinder#getNormal(primitives.Point3D)}.
	 */
	@Test
	public void testGetNormal() {
		 // ============ Equivalence Partitions Tests ==============
		  Cylinder t = new Cylinder(1,new Ray(new Point3D(1,0,0), new Vector(0,1,0)), 1);
	      assertEquals("Wrong normal to Cylinder", new Vector(1, 0, 0), t.getNormal(new Point3D(2, 0, 0)));
	}
	
	
	@Test
	public void testFindIntersections() 
	{
		
	}
	

}
