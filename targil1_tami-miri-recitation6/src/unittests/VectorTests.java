/**
 * 
 */
package unittests;

import static org.junit.Assert.*;

import org.junit.Test;
import static primitives.Util.*;

import primitives.Vector;

/**
 * Unit tests for primitives.Vector class
 */
public class VectorTests {

	/**
	 * Test method for {@link primitives.Vector#subtract(primitives.Vector)}.
	 */
	@Test
	public void testSubtract() {
		// ============ Equivalence Partitions Tests ==============
				 Vector v1 = new Vector(1, 1, 1);
			     Vector v2 = new Vector(2, 2, 2);
			     v1 = v1.subtract(v2);
			     assertEquals("Subtract() wrong result",(new Vector(-1,-1,-1)),v1);

			     Vector v3 = new Vector(3,3,3);
			     Vector v4 = new Vector(-1, -1, -1);
			     v3 = v3.subtract(v4);
			     assertEquals("Subtract() wrong result",new Vector(4,4,4),v3);
				
			 // =============== Boundary Values Tests ==================
		        try {
		            Vector v9=new Vector(1,1,1);
		            Vector v10=new Vector(1,1,1);
		            v9.subtract(v10);
		            fail("Vector (0,0,0) not valid");
		        }
		        catch  (IllegalArgumentException e)
		        {
		            assertTrue(e.getMessage()!= null);
		        }
	}

	/**
	 * Test method for {@link primitives.Vector#add(primitives.Vector)}.
	 */
	@Test
	public void testAdd() {
	 // ============ Equivalence Partitions Tests ==============
		 Vector v1 = new Vector(1, 1, 1);
	     Vector v2 = new Vector(1, 1, 1);
	     v1 = v1.add(v2);
	     assertEquals("add() wrong result",(new Vector(2, 2, 2)),v1);

	     Vector v3 = new Vector(-1, -1, -1);
	     Vector v4 = new Vector(-1, -1, -1);
	     v3 = v3.add(v4);
	     assertEquals("add() wrong result",new Vector(-2, -2, -2),v3);
		
	 // =============== Boundary Values Tests ==================
        try {
            Vector v9=new Vector(-1,-1,-1);
            Vector v10=new Vector(1,1,1);
            v9.add(v10);
            fail("Vector (0,0,0) not valid");
        }
        catch  (IllegalArgumentException e)
        {
            assertTrue(e.getMessage()!= null);
        }
	}

	/**
	 * Test method for {@link primitives.Vector#scale(double)}.
	 */
	@Test
	public void testScale() {
		 // ============ Equivalence Partitions Tests ==============
		 Vector v1 = new Vector(1, 1, 1);
	     v1 = v1.scale(3);
	     assertEquals("scale() wrong result",(new Vector(3, 3, 3)),v1);

	     Vector v3 = new Vector(-1, -1, -1);
	     v3 = v3.scale(-2);
	     assertEquals("scale() wrong result",new Vector(2,2,2),v3);
		
	 // =============== Boundary Values Tests ==================
	     try {
	            Vector v9=new Vector(-1,-1,-1);
	            v9.scale(0);
	            fail("Vector (0,0,0) not valid");
	        }
	        catch  (IllegalArgumentException e)
	        {
	            assertTrue(e.getMessage()!= null);
	        }
		
	}

	/**
	 * Test method for {@link primitives.Vector#dotProduct(primitives.Vector)}.
	 */
	@Test
	public void testDotProduct() {
		 // ============ Equivalence Partitions Tests ==============
		 Vector v1 = new Vector(1, 2, 3);
	     Vector v2 = new Vector(0,3,-2);
	     assertTrue("dotProduct() for orthogonal vectors is not zero",isZero(v1.dotProduct(v2))); 
	     Vector v3 = new Vector(-2, -4, -6);
	     assertTrue("dotProduct() wrong value",isZero(v1.dotProduct(v3) + 28));
	    
	}

	/**
	 * Test method for {@link primitives.Vector#crossProduct(primitives.Vector)}.
	 */
	@Test
	public void testCrossProduct() 
	{
        Vector v1 = new Vector(1, 2, 3);
        Vector v2 = new Vector(-2, -4, -6);

        // ============ Equivalence Partitions Tests ==============
        Vector v3 = new Vector(0, 3, -2);
        Vector vr = v1.crossProduct(v3);

        // Test that length of cross-product is proper (orthogonal vectors taken for simplicity)
        assertEquals("crossProduct() wrong result length", v1.length() * v3.length(), vr.length(), 0.00001);

        // Test cross-product result orthogonality to its operands
        assertTrue("crossProduct() result is not orthogonal to 1st operand", isZero(vr.dotProduct(v1)));
        assertTrue("crossProduct() result is not orthogonal to 2nd operand",isZero(vr.dotProduct(v3)));

        // =============== Boundary Values Tests ==================
        // test zero vector from cross-productof co-lined vectors
        try {
            v1.crossProduct(v2);
            fail("crossProduct() for parallel vectors does not throw an exception");
        } catch (Exception e) {}

	}

	/**
	 * Test method for {@link primitives.Vector#lengthSquared()}.
	 */
	@Test
	public void testLengthSquared() {
		// ============ Equivalence Partitions Tests ==============
		Vector v=new Vector(1,1,1);
		assertTrue("lengthSquard() doesn't return the right length",v.lengthSquared()==3);
		Vector v1=new Vector(-1,-1,-1);
		assertTrue("lengthSquard() doesn't return the right length",v1.lengthSquared()==3);
		
	}

	/**
	 * Test method for {@link primitives.Vector#length()}.
	 */
	@Test
	public void testLength() {
		// ============ Equivalence Partitions Tests ==============
		Vector v=new Vector(1,1,1);
		assertTrue("length() doesn't return the right length",v.length()==Math.sqrt(3));
		Vector v1=new Vector(-1,-1,-1);
		assertTrue("length() doesn't return the right length",v1.length()==Math.sqrt(3));
	}

	/**
	 * Test method for {@link primitives.Vector#normalize()}.
	 */
	@Test
	public void testNormalize() {
		 // ============ Equivalence Partitions Tests ==============
		Vector v=new Vector(3,0,0);
		v.normalize();
		assertEquals("normalize() result is not a unit vector",1,v.length(),1e-10);
		Vector v1=new Vector(-3,0,0);
		v1.normalize();
		assertEquals("normalize() result is not a unit vector",1,v1.length(),1e-10);
		
	}

	/**
	 * Test method for {@link primitives.Vector#normalized()}.
	 */
	@Test
	public void testNormalized() {
		// ============ Equivalence Partitions Tests ==============
		Vector v1=new Vector(3,0,0);
		Vector v2=v1.normalized();
		assertTrue("nurmalized() doesn't return new vector",!v2.equals(v1));
		assertEquals("normalized() result is not a unit vector",1,v2.length(),1e-10);
		Vector v3=new Vector(-3,0,0);
		v3=v3.normalized();
		assertEquals("normalized() result is not a unit vector",1,v3.length(),1e-10);
		 
				
	}

}
