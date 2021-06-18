// Copyright (C) 20080122 Taizo Kawano <tkawano at mshri.on.ca>
//
// This program is free software; you can redistribute it and/modify it 
// under the term of the GNU General Public License as published bythe Free Software Foundation;
// either version 2, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this file.  If not, write to the Free Software Foundation,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;
import ij.io.*;

public class Preprocessor_2 implements PlugInFilter {
	ImagePlus imp;

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		return DOES_ALL;
	}

	public void run(ImageProcessor ip) {
		ip.setInterpolate(true);
		boolean interpolate = PlotWindow.interpolate;
		Roi roi = imp.getRoi();
		
		/*---------------------------------------------------------------------*/
		// Rotation part. This is quoted from getIrregularProfile method within ProfilePlot,
		// and modified to get picture, not a line.
		/*---------------------------------------------------------------------*/

		//boolean calcXValues = cal!=null && cal.pixelWidth!=cal.pixelHeight;
		int n = ((PolygonRoi)roi).getNCoordinates();
		int[] x = ((PolygonRoi)roi).getXCoordinates();
		int[] y = ((PolygonRoi)roi).getYCoordinates();
		Rectangle r = roi.getBounds();
		int xbase = r.x;
		int ybase = r.y;
		double length = 0.0;
		double segmentLength;
		int xdelta, ydelta, iLength;
		double[] segmentLengths = new double[n];
		int[] dx = new int[n];
		int[] dy = new int[n];
		for (int i=0; i<(n-1); i++) {
			xdelta = x[i+1] - x[i];
			ydelta = y[i+1] - y[i];
			segmentLength = Math.sqrt(xdelta*xdelta+ydelta*ydelta);
			length += segmentLength;
			segmentLengths[i] = segmentLength;
			dx[i] = xdelta;
			dy[i] = ydelta;
		}
		double[] values = new double[(int)length];
		//make new ImageProcessor. newwidth = (int)length, newheight = 2*widthofstrip.
		//widthofstrip could be modified to fit practical image.
		int widthofstrip = 15;
		int newwidth = (int)length;
		int newheight = widthofstrip*2+1;
		ImageProcessor ip2 = new FloatProcessor(newwidth, newheight);

		//if (calcXValues) xValues = new float[(int)length];
		double leftOver = 1.0;
		double distance = 0.0;
		int index;
		double oldrx=0.0, oldry=0.0, xvalue=0.0;
		for (int i=0; i<n; i++) {
			double len = segmentLengths[i];
			if (len==0.0)
				continue;
			//To know the distance from point, need angle and theta angle of (vertical line.)
			double angle = Math.atan2(dy[i], dx[i]);
			double theta = angle-(Math.PI/2);
			//String anglestr = String.valueOf(angle*180/Math.PI);
			//String thetastr = String.valueOf(theta*180/Math.PI);
			//IJ.log(anglestr + " " + thetastr);
			double xinc = dx[i]/len;
			double yinc = dy[i]/len;
			double start = 1.0-leftOver;
			double rx = xbase+x[i]+start*xinc;
			double ry = ybase+y[i]+start*yinc;
			double len2 = len - start;
			int n2 = (int)len2;
			for (int j=0; j<=n2; j++) {
				index = (int)distance+j;
				//IJ.log(i+" "+index+" "+distance+" "+j);
				if (index<values.length) {
					//if (interpolate){
						//values[index] = ip.getInterpolatedValue(rx, ry);
						for (int k=0; k< newheight; k++){
							double drx = Math.cos(theta)* (widthofstrip-k);
							double dry = Math.sin(theta)* (widthofstrip-k);
							double rx2 = rx + drx;
							double ry2 = ry + dry;
							ip2.putPixelValue(index, k, ip.getInterpolatedValue(rx2, ry2));
							String xcood = String.valueOf(rx2);
							//IJ.log(xcood);
						}
					//}
					//else{
						//values[index] = ip.getPixelValue((int)(rx+0.5), (int)(ry+0.5));
						//IJ.log("not interpolate");
					//}
					/*if (calcXValues && index>0) {
						double deltax = cal.pixelWidth*(rx-oldrx);
						double deltay = cal.pixelHeight*(ry-oldry);
						xvalue += Math.sqrt(deltax*deltax + deltay*deltay);
						xValues[index]  = (float)xvalue;
					}*/
					oldrx = rx; oldry=ry;
				}
				rx += xinc;
				ry += yinc;
			}
			distance += len;
			leftOver = len2 - n2;
		}

		// return values;
		//ip2.resetMinAndMax();
		//here is just cropped image
		//ImageProcessor ip2short = ip2.convertToShort(false);
		//new ImagePlus("cropped", ip2short).show();
		//ImagePlus im2 = new ImagePlus("test", ip2);
		//new FileSaver(im2).saveAsText();
		
		/*---------------------------------------------------------------------*/
		//smoothing
		ImageProcessor ip3 = ip2.duplicate();
		ip3.medianFilter();
		ip3.smooth();
		//new ImagePlus("smooethen", ip3).show();

		/*---------------------------------------------------------------------*/
		//calculation of background: there should be better methods. but This is best with my ability for now.
		//You may just use subtract "background" of ImageJ default function.
		ImageProcessor ip4 = ip3.duplicate();
		int finallength = ip4.getWidth();
		float[][] floatarray = ip4.getFloatArray();
		//reference row number
		int refs = 5;
		for (int i=0; i < finallength; i++){
			//diff is linear interpolated difference of each pixel.
			float diff = (floatarray[i][refs-1]-floatarray[i][newheight-refs])/(newheight-refs*2+1);
			for(int j= refs; j < newheight-refs ;j++){				floatarray[i][j] = floatarray[i][j-1] - diff; 
				//smcolumn[j] = smcolumn[j-1]-diff;
			}
		}
		ip4.setFloatArray(floatarray);
		//ip4.resetMinAndMax();
		//new ImagePlus("background", ip4).show();
		
		/*---------------------------------------------------------------------*/
		//subtruct backgound from cropped image.
		ImageProcessor ip5 = ip2.duplicate();
		//ip5.copyBits(ip4, 0, 0, Blitter.DIFFERENCE);
		//first version uses difference. its not good because calculate absolute vale. use subtract.
		ip5.copyBits(ip4, 0, 0, Blitter.SUBTRACT);
		//ip5.resetMinAndMax();
		//convertToShort(false) clamped without scaling, but sets negative valuse to 0. this is different from octave subtracition method that I wrote before.
		ImageProcessor ip5short = ip5.convertToShort(false);
		ImagePlus im5 = new ImagePlus("subtracted", ip5short);
		//ImagePlus im5 = new ImagePlus("subtracted", ip5);
		im5.show();
		
		/*---------------------------------------------------------------------*/
		//our camera is 12bit, so save the 16 bit data.
		FileSaver fs = new FileSaver(im5);
		fs.saveAsText();
	}

}
