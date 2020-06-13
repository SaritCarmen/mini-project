package unittests;

import java.awt.Color;
import org.junit.Test;
import renderer.ImageWriter;


public class ImageWriterTest {

	    @Test
	    //create image in one color and grid in other color
	    public void writeToImage() {
	        ImageWriter imageWriter = new ImageWriter("new image", 1600, 1000, 800, 500);
	        int Nx = imageWriter.getNx();
	        int Ny = imageWriter.getNy();
	        for (int i = 0; i < Ny; i++) {
	            for (int j = 0; j < Nx; j++) {
	                if (i % 50 == 0 || j % 50 == 0) 
	                {
	                    imageWriter.writePixel(j, i, Color.RED);
	                } 
	                else 
	                {
	                    imageWriter.writePixel(j, i, Color.PINK);
	                }
	            }
	        }
	        imageWriter.writeToImage();
	    }

}
