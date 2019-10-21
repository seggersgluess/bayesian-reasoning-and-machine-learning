package Graphics;

import java.awt.Color;
import java.awt.Rectangle;

@SuppressWarnings("serial")
public class PieGraphics extends GraphicDevice{

	public static double [][] pieValues;
	public static Color [] pieColors;
	
	public void drawPie(Rectangle area) {
		
		int nValues = pieValues.length;
		
		if(nValues != pieColors.length) {
			System.out.println("Number of values and colors are not identical.");
			return;
		}
		
		double total = 0.0D;
		for (int i = 0; i < nValues; i++) {
			total += pieValues[i][0];
		}
		double curValue = 0.0D;
		int startAngle = 0;
		for (int i = 0; i < nValues; i++) {
			startAngle = (int) (curValue * 360 / total);
		    int arcAngle = (int) (pieValues[i][0] * 360 / total);
		    g2.setColor(pieColors[i]);
		    g2.fillArc(area.x, area.y, area.width, area.height, 
		    startAngle, arcAngle);
		    curValue += pieValues[i][0];
		}
	}

	
	public static void setPieValues(double [][] values) {
		pieValues = values;
	}
	
	
	public static void setPieColors(Color [] colors) {
		pieColors = colors;
	}
	
}
