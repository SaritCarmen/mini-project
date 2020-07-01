package unittests;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import elements.AmbientLight;
import elements.Camera;
import elements.PointLight;
import elements.SpotLight;
import geometries.Cylinder;
import geometries.Geometries;
import geometries.Plane;
import geometries.Sphere;

import k_means.Cluster;
import k_means.Kmeans;
import k_means.Point;
import primitives.Color;
import primitives.Material;
import primitives.Point3D;
import primitives.Ray;
import primitives.Vector;
import renderer.ImageWriter;
import renderer.Render;
import scene.Scene;


public class miniProject2 {

	/**
	 * implement the mini project2 part2-use Geometries objects to create the hierarchy
	 */
	@Test
	public void miniPtest1() {
		Scene scene = new Scene("Test scene");
        scene.setCamera(new Camera(new Point3D(0, 0, -1000), new Vector(0, 0, 1), new Vector(0, -1, 0)));
        scene.setDistance(1000);
        scene.setBackground(new Color(java.awt.Color.BLACK));
        scene.setAmbientLight(new AmbientLight(new Color(java.awt.Color.WHITE), 0.15));
        
        Sphere A=new Sphere(new Color(255,0,0), new Material(0.001, 1.5, 100, 0, 0.7), 20, new Point3D(0, 0, 0));
        Sphere B=new Sphere(new Color(100,0,50), new Material(0.001, 1.5, 100, 0, 0.7), 20, new Point3D(60, 0, 400));
        Sphere C=new Sphere(new Color(56,115,8), new Material(0.001, 1.5, 100, 0, 0.7), 20, new Point3D(30,0,200));
        
       
        Sphere G=  new Sphere(new Color(100,0,100), new Material(0.001, 1.5, 100, 0.2, 0.7), 5, new Point3D(-60, 15, 0));
        Sphere H=  new Sphere(new Color(70,50,0), new Material(0.001, 1.5, 100, 0, 1), 25, new Point3D(-70,-5,80));
        Sphere I=  new Sphere(new Color(java.awt.Color.BLUE), new Material(0.001, 1.5, 100, 0.2, 0.7), 10, new Point3D(-40, 10, 0));
        Sphere J=  new Sphere(new Color(100,0,100), new Material(0.001, 1.5, 100, 0.2, 0.7), 5, new Point3D(70, 15, 0));
        
        
        Plane P=  new Plane(new Color(java.awt.Color.BLACK), new Material(0.5, 0.5, 0, 0.5, 0.7), new Point3D(0, 20, 0), new Point3D(4, 20, 0),new Point3D(4,20, 2));
        
        
        Geometries geos2=new Geometries(A,B,C);
        Geometries geos3=new Geometries(J);
        Geometries geos4= new Geometries(G,H,I);
        Geometries geos1=new Geometries(geos2,geos3,geos4,P);
               scene.addGeometries(geos1);
        
      
      scene.addLights(new SpotLight(new Color(700, 400, 400),new Vector(-1, 1, 4),new Point3D(40, -40, -115),  1, 4E-4, 2E-5,10));
      
      scene.addLights(
   		   new PointLight(new Color(700, 400, 400), new Point3D(0, -10,0),1, 0.001, 0.0005,10));
   
     scene.addLights(new SpotLight(new Color(700, 400, 400), //
				new Vector(0, 0, 1), new Point3D(60, -50, 0), 1, 4E-5, 2E-7,10));
           
                    
        ImageWriter imageWriter = new ImageWriter("miniP2 part1", 200, 200, 500, 500);
        Render render = new Render(imageWriter, scene).setMultithreading(3).setDebugPrint();//60's without multithreading and 45's with.

        
        render.renderImage();
        render.writeToImage();
	}
	
	/**
	 * implement the mini project2 part3- using an algorithm (k means) that automatically re-arranges the “flat” scene’s geometries set into geometries hierarchy 
	 */
	@Test
	public void miniPtest2() {
		Scene scene = new Scene("Test scene");
        scene.setCamera(new Camera(new Point3D(0, 0, -1000), new Vector(0, 0, 1), new Vector(0, -1, 0)));
        scene.setDistance(1000);
        scene.setBackground(Color.BLACK);
        scene.setAmbientLight(new AmbientLight(new Color(java.awt.Color.WHITE), 0.15));
        
        
        Sphere A=new Sphere(new Color(255,0,0), new Material(0.001, 1.5, 100, 0, 0.7), 20, new Point3D(0, 0, 0));
        Sphere B=new Sphere(new Color(100,0,50), new Material(0.001, 1.5, 100, 0, 0.7), 20, new Point3D(60, 0, 400));
        Sphere C=new Sphere(new Color(56,115,8), new Material(0.001, 1.5, 100, 0, 0.7), 20, new Point3D(30,0,200));
        
       
        Sphere G=  new Sphere(new Color(100,0,100), new Material(0.001, 1.5, 100, 0.2, 0.7), 5, new Point3D(-60, 15, 0));
        Sphere H=  new Sphere(new Color(70,50,0), new Material(0.001, 1.5, 100, 0, 1), 25, new Point3D(-70,-5,80));
        Sphere I=  new Sphere(new Color(java.awt.Color.BLUE), new Material(0.001, 1.5, 100, 0.2, 0.7), 10, new Point3D(-40, 10, 0));
        Sphere J=  new Sphere(new Color(100,0,100), new Material(0.001, 1.5, 100, 0.2, 0.7), 5, new Point3D(70, 15, 0));
        
        
        Plane P=  new Plane(new Color(java.awt.Color.BLACK), new Material(0.5, 0.5, 0, 0.5, 0.7), new Point3D(0, 20, 0), new Point3D(4, 20, 0),new Point3D(4,20, 2));
        
        List<Point> points=new ArrayList<Point>();
		points.add(new Point(A));
		points.add(new Point(B));
		points.add(new Point(C));		
		points.add(new Point(J));
		points.add(new Point(G));
		points.add(new Point(H));
		points.add(new Point(I));
		points.add(new Point(P));
		
		Kmeans kmeans = new Kmeans();
        kmeans.init(points);
        kmeans.calculate();
        List<Cluster> clusters=kmeans.getClusters();
        for(Cluster c: clusters) {
        	Geometries geos = new Geometries();       	
        	for(Point p: c.getPoints()) {
        		geos.add(p.getGeometry());
        	}
        	scene.addGeometries(geos);
        }
      
      scene.addLights(new SpotLight(new Color(700, 400, 400),new Vector(-1, 1, 4),new Point3D(40, -40, -115),  1, 4E-4, 2E-5,10));
      
      scene.addLights(
   		   new PointLight(new Color(700, 400, 400), new Point3D(0, -10,0),1, 0.001, 0.0005,10));
   
    scene.addLights(new SpotLight(new Color(700, 400, 400), //
				new Vector(0, 0, 1), new Point3D(60, -50, 0), 1, 4E-5, 2E-7,10));
           
                    
        ImageWriter imageWriter = new ImageWriter("miniP2 part2", 200, 200, 500, 500);
        Render render = new Render(imageWriter, scene).setMultithreading(3).setDebugPrint();//60's without multithreading and 45's with.

        
        render.renderImage();
        render.writeToImage();
	}


}