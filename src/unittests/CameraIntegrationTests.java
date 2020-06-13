package unittests;

import static org.junit.Assert.*;

import java.util.List;


import org.junit.Test;

import elements.Camera;
import geometries.Intersectable.GeoPoint;
import geometries.Plane;
import geometries.Sphere;
import geometries.Triangle;
import primitives.Point3D;
import primitives.Vector;

public class CameraIntegrationTests {
	Camera camera1=new Camera(Point3D.ZERO,new Vector(0,0,1),new Vector(0,-1,0));
	Camera camera2=new Camera(new Point3D(0,0,-0.5),new Vector(0,0,1),new Vector(0,-1,0));
	
	/**
	 * test sphere is against the view plane(2 points)
	 */
	@Test
	public void constructRayThroughPixelWithSphere1()
	{
		Sphere sphere =  new Sphere(1, new Point3D(0, 0, 3));
        List<GeoPoint> results;
        int count = 0;
        int Nx =3;
        int Ny =3;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                results = sphere.findIntersections(camera1.constructRayThroughPixel(Nx, Ny, j, i, 1, 3, 3));
                if (results != null)
                    count += results.size();
            }
        }

        assertEquals("Bad ray",2,count);
        System.out.println("count: "+count);

		
	}
	/**
	 * the shpere starts befor the view plane (18 points)
	 */
	@Test
	public void constructRayThroughPixelWithSphere2() {
		Sphere sphere =  new Sphere(2.5, new Point3D(0, 0, 2.5));
        List<GeoPoint> results;
        int count = 0;
        int Nx =3;
        int Ny =3;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                results = sphere.findIntersections(camera2.constructRayThroughPixel(Nx, Ny, j, i, 1, 3, 3));
                if (results != null)
                    count += results.size();
            }
        }

        assertEquals("Bad ray",18,count);
        System.out.println("count: "+count);

		
	}
	/**
	 * five rayes cross the sphere (10 points)
	 */
	@Test
	public void constructRayThroughPixelWithSphere3() {
		Sphere sphere =  new Sphere(2, new Point3D(0, 0, 2));
        List<GeoPoint> results;
        int count = 0;
        int Nx =3;
        int Ny =3;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                results = sphere.findIntersections(camera2.constructRayThroughPixel(Nx, Ny, j, i, 1, 3, 3));
                if (results != null)
                    count += results.size();
            }
        }

        assertEquals("Bad ray",10,count);
        System.out.println("count: "+count);

	}
	
	/**
	 * the camera is in the shpere (9 points)
	 */
	@Test
	public void constructRayThroughPixelWithSphere4() {
		Sphere sphere =  new Sphere(4, new Point3D(0, 0, 0));
		List<GeoPoint> results;
        int count = 0;
        int Nx =3;
        int Ny =3;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                results = sphere.findIntersections(camera1.constructRayThroughPixel(Nx, Ny, j, i, 1, 3, 3));
                if (results != null)
                    count += results.size();
            }
        }

        assertEquals("Bad ray",9,count);
        System.out.println("count: "+count);

	}
	
	/**
	 * the sphere is befor the camera(0 points)
	 */
	@Test
	public void constructRayThroughPixelWithSphere5() {
		Sphere sphere =  new Sphere(0.5, new Point3D(0, 0, -1));
		List<GeoPoint> results;
        int count = 0;
        int Nx =3;
        int Ny =3;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                results = sphere.findIntersections(camera1.constructRayThroughPixel(Nx, Ny, j, i, 1, 3, 3));
                if (results != null)
                    count += results.size();
            }
        }

        assertEquals("Bad ray",0,count);
        System.out.println("count: "+count);

	}
	
	/**
	 * the plane is orthogonal to the view plane (9 points)
	 */
	@Test
    public void constructRayThroughPixelWithPlane1() {
		Plane plane =  new Plane( new Point3D(0, 0, 2), new Vector(0,0,1));
		List<GeoPoint> results;
        int count = 0;
        int Nx =3;
        int Ny =3;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                results = plane.findIntersections(camera1.constructRayThroughPixel(Nx, Ny, j, i, 1, 3, 3));
                if (results != null)
                    count += results.size();
            }
        }

        assertEquals("Bad ray",9,count);
        System.out.println("count: "+count);

	}
	
	/**
	 * the plane is at an angle to the view plane (9 points)
	 */
	@Test
    public void constructRayThroughPixelWithPlane2() {
		Plane plane =  new Plane( new Point3D(0, 0, 2), new Point3D(1,-3,1), new Point3D(0,-6,0));
		List<GeoPoint> results;
        int count = 0;
        int Nx =3;
        int Ny =3;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                results = plane.findIntersections(camera1.constructRayThroughPixel(Nx, Ny, j, i, 1, 3, 3));
                if (results != null)
                    count += results.size();
            }
        }

        assertEquals("Bad ray",9,count);
        System.out.println("count: "+count);
	}
	
	/**
	 *  the plane is at an angle to the view plane but cross the top of the view plane (6 points)
	 */
	@Test
    public void constructRayThroughPixelWithPlane3() {
		Plane plane =  new Plane( new Point3D(0, 0, 5), new Point3D(0,-5,0), new Point3D(1,-2.5,2.5));
		List<GeoPoint> results;
        int count = 0;
        int Nx =3;
        int Ny =3;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                results = plane.findIntersections(camera1.constructRayThroughPixel(Nx, Ny, j, i, 1, 3, 3));
                if (results != null)
                    count += results.size();
            }
        }

        assertEquals("Bad ray",6,count);
        System.out.println("count: "+count);
	}
	
	/**
	 * the triangle against only one of the pixels
	 */
	@Test
    public void constructRayThroughPixelWithTriangle1() {
		Triangle triangle =  new Triangle( new Point3D(0, -1, 2), new Point3D(1,1,2), new Point3D(-1,1,2));
		List<GeoPoint> results;
        int count = 0;
        int Nx =3;
        int Ny =3;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                results = triangle.findIntersections(camera1.constructRayThroughPixel(Nx, Ny, j, i, 1, 3, 3));
                if (results != null)
                    count += results.size();
            }
        }

        assertEquals("Bad ray",1,count);
        System.out.println("count: "+count);
	}
	
	/**
	 * the triangle against two of the pixels
	 */
	@Test
    public void constructRayThroughPixelWithTriangle2() {
		Triangle triangle =  new Triangle( new Point3D(0, -20, 2), new Point3D(1,1,2), new Point3D(-1,1,2));
		List<GeoPoint> results;
        int count = 0;
        int Nx =3;
        int Ny =3;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                results = triangle.findIntersections(camera1.constructRayThroughPixel(Nx, Ny, j, i, 1, 3, 3));
                if (results != null)
                    count += results.size();
            }
        }

        assertEquals("Bad ray",2,count);
        System.out.println("count: "+count);
	}
}
